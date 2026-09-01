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
    target_concentration_m: float
    thf_count: int
    amine_count_initial: int
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
            "li_electron_pairs": self.li_electron_pairs,
            "temperature_k": self.temperature_k,
            "pressure_bar": self.pressure_bar,
        }


def make_system_spec(
    system_id: str, campaign: dict[str, Any], systems: dict[str, Any]
) -> SystemSpec:
    amine, target = parse_system_id(system_id)
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
            target_concentration_m=target,
            thf_count=thf_count,
            thf_molar_volume_l_mol=float(systems["thf"]["molar_volume_l_mol"]),
            amine_molar_volume_l_mol=float(amines[amine]["molar_volume_l_mol"]),
        )
    return SystemSpec(
        system_id=system_id,
        amine=amine,
        target_concentration_m=target,
        thf_count=thf_count,
        amine_count_initial=amine_count,
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
    thresholds = methods.get("localization_thresholds", {})
    if float(thresholds.get("electron_count_min", 2)) >= float(
        thresholds.get("electron_count_max", 0)
    ):
        errors.append("electron_count_min must be smaller than electron_count_max")
    return errors
