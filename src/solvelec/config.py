from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .composition import initial_amine_count


class ConfigurationError(ValueError):
    """Raised when a campaign or method configuration is unsafe or malformed."""


def repository_root(start: Path | None = None) -> Path:
    current = (start or Path.cwd()).resolve()
    for candidate in (current, *current.parents):
        if (candidate / "configs" / "campaign.yaml").is_file():
            return candidate
    raise ConfigurationError("Could not locate repository root containing configs/campaign.yaml")


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a config written in the JSON subset of YAML.

    Keeping committed configs in this subset lets the core CLI validate and
    inspect campaigns before PyYAML or Snakemake is installed.
    """

    config_path = Path(path)
    try:
        with config_path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except json.JSONDecodeError as exc:
        raise ConfigurationError(
            f"{config_path} must remain valid JSON (which is also valid YAML): {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise ConfigurationError(f"{config_path} must contain an object at the top level")
    return value


def load_repository_configs(root: Path | None = None) -> tuple[dict[str, Any], ...]:
    base = root or repository_root()
    return (
        load_config(base / "configs" / "campaign.yaml"),
        load_config(base / "configs" / "systems.yaml"),
        load_config(base / "configs" / "methods.yaml"),
    )


_SYSTEM_RE = re.compile(r"^(?P<amine>[a-z0-9]+)_(?P<value>[0-9]+(?:p[0-9]+)?)m$")


def parse_system_id(system_id: str) -> tuple[str | None, float]:
    if system_id == "pure_thf":
        return None, 0.0
    match = _SYSTEM_RE.fullmatch(system_id)
    if not match:
        raise ConfigurationError(
            f"Invalid system id {system_id!r}; expected pure_thf or <amine>_<value>m"
        )
    return match.group("amine"), float(match.group("value").replace("p", "."))


def format_system_id(amine: str | None, concentration_m: float) -> str:
    if amine is None or concentration_m == 0:
        return "pure_thf"
    value = f"{concentration_m:g}".replace(".", "p")
    return f"{amine}_{value}m"


@dataclass(frozen=True)
class SystemSpec:
    system_id: str
    amine: str | None
    target_concentration_m: float | None
    thf_count: int
    amine_count_initial: int
    component_counts: dict[str, int]
    composition_basis: str
    target_amine_mole_fraction: float
    li_electron_pairs: int
    temperature_k: float
    pressure_bar: float

    def as_dict(self) -> dict[str, Any]:
        return {
            "system_id": self.system_id,
            "amine": self.amine,
            "target_concentration_m": self.target_concentration_m,
            "thf_count": self.thf_count,
            "amine_count_initial": self.amine_count_initial,
            "component_counts": dict(self.component_counts),
            "composition_basis": self.composition_basis,
            "target_amine_mole_fraction": self.target_amine_mole_fraction,
            "li_electron_pairs": self.li_electron_pairs,
            "temperature_k": self.temperature_k,
            "pressure_bar": self.pressure_bar,
        }


def make_system_spec(
    system_id: str, campaign: dict[str, Any], systems: dict[str, Any]
) -> SystemSpec:
    explicit = systems.get("compositions", {}).get(system_id)
    if explicit is not None:
        components = explicit.get("components", {})
        if not isinstance(components, dict) or not components:
            raise ConfigurationError(
                f"Explicit composition {system_id!r} must define component counts"
            )
        component_counts = {str(name): int(count) for name, count in components.items()}
        known = {"thf", *systems.get("amines", {})}
        unknown = sorted(set(component_counts) - known)
        if unknown:
            raise ConfigurationError(
                f"Explicit composition {system_id!r} has unknown components: {unknown}"
            )
        if any(count <= 0 for count in component_counts.values()):
            raise ConfigurationError(f"Explicit composition {system_id!r} counts must be positive")
        amines = sorted(set(component_counts) - {"thf"})
        if len(amines) > 1:
            raise ConfigurationError(
                f"Explicit composition {system_id!r} currently supports one amine"
            )
        amine = amines[0] if amines else None
        thf_count = component_counts.get("thf", 0)
        amine_count = component_counts.get(amine, 0) if amine else 0
        target: float | None = None
        composition_basis = str(explicit.get("composition_basis", "explicit_molecule_counts"))
        total_count = sum(component_counts.values())
        target_mole_fraction = float(
            explicit.get("target_amine_mole_fraction", amine_count / total_count)
        )
        if abs(target_mole_fraction - amine_count / total_count) > 1.0e-12:
            raise ConfigurationError(
                f"Explicit composition {system_id!r} target mole fraction disagrees with counts"
            )
    else:
        amine, parsed_target = parse_system_id(system_id)
        target = parsed_target
        thf_count = int(campaign["thf_count"])
        if thf_count <= 0:
            raise ConfigurationError("thf_count must be positive")
        if amine is None:
            amine_count = 0
        else:
            amines = systems.get("amines", {})
            if amine not in amines:
                raise ConfigurationError(f"Unknown amine {amine!r} in system {system_id!r}")
            amine_count = initial_amine_count(
                target_concentration_m=parsed_target,
                thf_count=thf_count,
                thf_molar_volume_l_mol=float(systems["thf"]["molar_volume_l_mol"]),
                amine_molar_volume_l_mol=float(amines[amine]["molar_volume_l_mol"]),
            )
        component_counts = {"thf": thf_count}
        if amine:
            component_counts[amine] = amine_count
        composition_basis = "target_molarity"
        target_mole_fraction = amine_count / sum(component_counts.values())
    return SystemSpec(
        system_id=system_id,
        amine=amine,
        target_concentration_m=target,
        thf_count=thf_count,
        amine_count_initial=amine_count,
        component_counts=component_counts,
        composition_basis=composition_basis,
        target_amine_mole_fraction=target_mole_fraction,
        li_electron_pairs=int(campaign["li_electron_pairs"]),
        temperature_k=float(campaign["temperature_k"]),
        pressure_bar=float(campaign["pressure_bar"]),
    )


def campaign_matrix(
    campaign_name: str, campaign: dict[str, Any], systems: dict[str, Any]
) -> list[tuple[SystemSpec, int]]:
    campaigns = campaign.get("campaigns", {})
    if campaign_name not in campaigns:
        raise ConfigurationError(
            f"Unknown campaign {campaign_name!r}; choose from {sorted(campaigns)}"
        )
    definition = campaigns[campaign_name]
    replicas = definition.get("replicas", campaign.get("replicas", [1]))
    rows: list[tuple[SystemSpec, int]] = []
    for system_id in definition["systems"]:
        spec = make_system_spec(system_id, campaign, systems)
        for replica in replicas:
            if int(replica) <= 0:
                raise ConfigurationError("Replica identifiers must be positive integers")
            rows.append((spec, int(replica)))
    return rows


def validate_repository_configs(root: Path | None = None) -> list[str]:
    base = root or repository_root()
    campaign, systems, methods = load_repository_configs(base)
    errors: list[str] = []

    if campaign.get("schema_version") != 1:
        errors.append("campaign.schema_version must be 1")
    for name in campaign.get("campaigns", {}):
        try:
            campaign_matrix(name, campaign, systems)
        except (ConfigurationError, KeyError, TypeError, ValueError) as exc:
            errors.append(f"campaign {name}: {exc}")
    if not methods.get("cp2k", {}).get("uks", False):
        errors.append("methods.cp2k.uks must be true for the one-electron doublet")
    if methods.get("cp2k", {}).get("production_status") != (
        "PROVISIONAL_UNTIL_BASIS_FUNCTIONAL_GATE"
    ):
        errors.append("CP2K production status must remain explicitly provisional until G4")
    classical = methods.get("classical_md", {})
    for key in (
        "timestep_fs",
        "nvt_ns",
        "npt_equilibration_ns",
        "production_ns",
        "checkpoint_minutes",
        "trajectory_stride_ps",
    ):
        try:
            if float(classical[key]) <= 0:
                raise ValueError
        except (KeyError, TypeError, ValueError):
            errors.append(f"methods.classical_md.{key} must be positive")
    classical_validation = methods.get("classical_validation", {})
    for key in (
        "minimum_trajectory_fraction",
        "density_half_relative_tolerance",
        "replica_density_relative_tolerance",
    ):
        try:
            value = float(classical_validation[key])
            if not 0 < value <= 1:
                raise ValueError
        except (KeyError, TypeError, ValueError):
            errors.append(f"methods.classical_validation.{key} must be in (0, 1]")
    trajectory_analysis = methods.get("trajectory_analysis", {})
    for key in (
        "analysis_stride_ps",
        "minimum_effective_samples",
        "rdf_bin_width_angstrom",
        "rdf_max_angstrom",
        "eda_thf_contact_cutoff_angstrom",
        "hydrogen_bond_distance_angstrom",
        "cavity_hbond_shell_thickness_angstrom",
        "cavity_visualization_radius_angstrom",
        "minimum_snapshot_separation_ps",
        "decorrelation_multiplier",
    ):
        try:
            if float(trajectory_analysis[key]) <= 0:
                raise ValueError
        except (KeyError, TypeError, ValueError):
            errors.append(f"methods.trajectory_analysis.{key} must be positive")
    for key in ("minimum_analysis_frames", "void_grid_points_per_axis"):
        try:
            value = int(trajectory_analysis[key])
            if value < 2 or value != float(trajectory_analysis[key]):
                raise ValueError
        except (KeyError, TypeError, ValueError):
            errors.append(f"methods.trajectory_analysis.{key} must be an integer >= 2")
    try:
        refinement_levels = int(trajectory_analysis["void_refinement_levels"])
        if refinement_levels < 0 or refinement_levels != float(
            trajectory_analysis["void_refinement_levels"]
        ):
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append("methods.trajectory_analysis.void_refinement_levels must be an integer >= 0")
    try:
        if int(trajectory_analysis["snapshots_per_replica"]) != 1:
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append("methods.trajectory_analysis.snapshots_per_replica must be 1 for stage A")
    try:
        equilibrated_start = float(trajectory_analysis["equilibrated_start_ns"])
        if not 0 <= equilibrated_start < float(classical["production_ns"]):
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append(
            "methods.trajectory_analysis.equilibrated_start_ns must be in [0, production_ns)"
        )
    try:
        angle = float(trajectory_analysis["hydrogen_bond_angle_degree"])
        if not 0 < angle <= 180:
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append("methods.trajectory_analysis.hydrogen_bond_angle_degree must be in (0, 180]")
    try:
        if float(trajectory_analysis["cavity_hbond_bridge_margin_angstrom"]) < 0:
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append(
            "methods.trajectory_analysis.cavity_hbond_bridge_margin_angstrom must be non-negative"
        )
    stage_b = methods.get("stage_b", {})
    for key in (
        "candidate_site_count",
        "void_grid_points_per_axis",
        "smoke_replicas_per_system",
    ):
        try:
            value = int(stage_b[key])
            if value < 1 or value != float(stage_b[key]):
                raise ValueError
        except (KeyError, TypeError, ValueError):
            errors.append(f"methods.stage_b.{key} must be a positive integer")
    try:
        refinement = int(stage_b["void_refinement_levels"])
        if refinement < 0 or refinement != float(stage_b["void_refinement_levels"]):
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append("methods.stage_b.void_refinement_levels must be an integer >= 0")
    for key in ("minimum_site_separation_angstrom", "minimum_surface_clearance_angstrom"):
        try:
            value = float(stage_b[key])
            if value <= 0:
                raise ValueError
        except (KeyError, TypeError, ValueError):
            errors.append(f"methods.stage_b.{key} must be positive")
    candidate_ids: list[str] = []
    try:
        definitions = stage_b["candidate_pairs"]
        if not isinstance(definitions, list) or not definitions:
            raise ValueError
        for definition in definitions:
            candidate_id = str(definition["id"])
            if not re.fullmatch(r"[a-z][a-z0-9_-]*", candidate_id):
                raise ValueError
            if float(definition["target_li_cavity_distance_angstrom"]) <= 0:
                raise ValueError
            if float(definition["tolerance_angstrom"]) <= 0:
                raise ValueError
            candidate_ids.append(candidate_id)
        if len(candidate_ids) != len(set(candidate_ids)):
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append("methods.stage_b.candidate_pairs must define unique positive candidates")
    if str(stage_b.get("smoke_candidate_id", "")) not in candidate_ids:
        errors.append("methods.stage_b.smoke_candidate_id must name a configured candidate")
    stage_b_smoke = methods.get("stage_b_smoke", {})
    if stage_b_smoke.get("scientific_status") != "NUMERICAL_SMOKE_ONLY":
        errors.append("Stage-B CP2K smoke status must remain explicitly non-production")
    try:
        smoke_li_valence = int(stage_b_smoke["li_pseudopotential_valence_electrons"])
        smoke_li_target = float(stage_b_smoke["li_target_valence_electrons"])
        if smoke_li_valence != 3 or smoke_li_target != 2.0:
            raise ValueError
        if stage_b_smoke["li_potential"] != "GTH-PBE-q3":
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append("Stage-B Li+ must pair GTH-PBE-q3 with a 2.0-electron cDFT target")
    try:
        mechanism_states = stage_b_smoke["mechanism_states"]
        if not isinstance(mechanism_states, list) or len(mechanism_states) != 2:
            raise ValueError
        state_targets = {
            str(record["id"]): float(record["li_target_valence_electrons"])
            for record in mechanism_states
        }
        if state_targets != {"li0_diabatic": 3.0, "li_plus_e_diabatic": 2.0}:
            raise ValueError
        if any(not str(record.get("interpretation", "")).strip() for record in mechanism_states):
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append("Stage-B mechanism smoke must define interpreted Li0=3.0 and Li+e=2.0 states")
    for key in ("cutoff_ry", "rel_cutoff_ry", "max_scf", "cube_stride"):
        try:
            if float(stage_b_smoke[key]) <= 0:
                raise ValueError
        except (KeyError, TypeError, ValueError):
            errors.append(f"methods.stage_b_smoke.{key} must be positive")
    try:
        cdft_eps_scf = float(stage_b_smoke["cdft_eps_scf"])
        if not 0 < cdft_eps_scf <= 5.0e-2:
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append("methods.stage_b_smoke.cdft_eps_scf must be in (0, 0.05]")
    cp2k = methods.get("cp2k", {})
    try:
        if cp2k["li_potential"] != "GTH-PBE-q3":
            raise ValueError
        if int(cp2k["li_pseudopotential_valence_electrons"]) != 3:
            raise ValueError
        if float(cp2k["li_target_valence_electrons"]) != 2.0:
            raise ValueError
    except (KeyError, TypeError, ValueError):
        errors.append("production CP2K Li+ must pair GTH-PBE-q3 with target 2.0")
    thresholds = methods.get("localization_thresholds", {})
    if float(thresholds.get("electron_count_min", 2)) >= float(
        thresholds.get("electron_count_max", 0)
    ):
        errors.append("electron_count_min must be smaller than electron_count_max")
    return errors
