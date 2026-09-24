#!/usr/bin/env python3
"""Stage-B3 N-H association diagnostics and a minimal amine-series pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from solvelec.candidates import read_xyz, write_xyz
from solvelec.cube import read_cube
from solvelec.localization_analysis import (
    candidate_cube_problems,
    cube_compatibility_problems,
    infer_molecular_topology,
    read_cp2k_cell,
    state_localization_record,
)
from solvelec.nh_environment import (
    infer_nh_donors,
    nh_environment,
    select_nh_contrast_sites,
)
from solvelec.parsers import HARTREE_TO_EV, parse_cp2k_text
from solvelec.preferential_solvation import local_composition, molecular_heavy_atom_centers
from solvelec.provenance import sha256_file
from solvelec.rendering import render_stage_b3_amine_series_cp2k

SCIENTIFIC_STATUS = "NUMERICAL_AMINE_SERIES_NH_ASSOCIATION_SCREEN_ONLY"


def _read_json(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def _write_json(path: str | Path, value: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)


def _component_counts_from_spec(spec: dict[str, Any]) -> dict[str, int]:
    configured = spec.get("component_counts")
    if isinstance(configured, dict) and configured:
        counts = {str(name): int(value) for name, value in configured.items()}
    else:
        counts: dict[str, int] = {}
        if int(spec.get("thf_count", 0)):
            counts["thf"] = int(spec["thf_count"])
        if spec.get("amine") and int(spec.get("amine_count_initial", 0)):
            counts[str(spec["amine"])] = int(spec["amine_count_initial"])
    if not counts or any(value <= 0 for value in counts.values()):
        raise ValueError("Stage-B3 spec has no positive solvent component count")
    return counts


def _component_definitions(systems: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        "thf": systems["thf"],
        **{str(name): value for name, value in systems.get("amines", {}).items()},
    }


def _validate_checksum_gate(summary_path: str | Path, gate_path: str | Path) -> None:
    summary = Path(summary_path)
    gate = Path(gate_path)
    if not summary.is_file() or not gate.is_file():
        raise ValueError(f"missing immutable Stage-B3 handoff: {summary} / {gate}")
    if not _read_json(summary).get("ready"):
        raise ValueError(f"handoff summary is not ready: {summary}")
    digest = hashlib.sha256(summary.read_bytes()).hexdigest()
    if gate.read_text(encoding="utf-8").split() != ["sha256", digest, summary.name]:
        raise ValueError(f"handoff checksum gate does not match {summary}")


def _source_structure(
    manifest: dict[str, Any],
) -> tuple[list[str], np.ndarray, Path, np.ndarray, dict[str, Any]]:
    candidates = [record for record in manifest.get("candidates", []) if record.get("ready")]
    if not manifest.get("ready") or not candidates:
        raise ValueError("Stage-B3 source candidate manifest is not ready")
    structure = candidates[0]["structure"]
    xyz_path = Path(structure["xyz"]["path"])
    cell_path = Path(structure["cell"]["path"])
    if sha256_file(xyz_path) != structure["xyz"]["sha256"]:
        raise ValueError("Stage-B3 source candidate XYZ checksum mismatch")
    if sha256_file(cell_path) != structure["cell"]["sha256"]:
        raise ValueError("Stage-B3 source candidate cell checksum mismatch")
    elements, positions, _comment = read_xyz(xyz_path)
    solvent_indices = [
        index for index, element in enumerate(elements) if element.upper() not in {"LI", "GH"}
    ]
    if len(solvent_indices) + 2 != len(elements):
        raise ValueError("Stage-B3 source must contain exactly one Li and one ghost center")
    solvent_elements = [elements[index] for index in solvent_indices]
    solvent_positions = np.asarray(positions[solvent_indices], dtype=float)
    return solvent_elements, solvent_positions, cell_path, read_cp2k_cell(cell_path), structure


def _topology_and_donors(
    elements: list[str],
    positions: np.ndarray,
    cell: np.ndarray,
    spec: dict[str, Any],
    systems: dict[str, Any],
    methods: dict[str, Any],
):
    settings = methods["stage_b3_amine_series"]
    topology = infer_molecular_topology(
        elements,
        positions,
        cell,
        _component_counts_from_spec(spec),
        _component_definitions(systems),
        bond_scale=float(methods["stage_b_localization"]["covalent_bond_scale"]),
    )
    if topology.problems:
        raise ValueError("Stage-B3 molecular topology failed: " + "; ".join(topology.problems))
    donors = infer_nh_donors(
        elements,
        positions,
        cell,
        topology,
        maximum_nh_bond_angstrom=float(settings["maximum_nh_bond_angstrom"]),
    )
    return topology, donors


def _environment(
    point: list[float] | np.ndarray,
    donors,
    positions: np.ndarray,
    cell: np.ndarray,
    settings: dict[str, Any],
) -> dict[str, Any]:
    return nh_environment(
        point,
        donors,
        positions,
        cell,
        shell_radii_angstrom=settings["nh_shell_radii_angstrom"],
        inward_angle_minimum_degrees=float(settings["inward_angle_minimum_degrees"]),
        distance_decay_angstrom=float(settings["distance_decay_angstrom"]),
    )


def run_reanalyze(args: argparse.Namespace) -> int:
    analysis = _read_json(args.analysis)
    methods = _read_json(args.methods)
    systems = _read_json(args.systems)
    settings = methods["stage_b3_amine_series"]
    if not analysis.get("ready"):
        raise ValueError("existing Stage-B2C analysis is not ready")
    seed_input = analysis.get("inputs", {}).get("seed_metadata", {})
    metadata_path = Path(seed_input.get("path", ""))
    if not metadata_path.is_file() or sha256_file(metadata_path) != seed_input.get("sha256"):
        raise ValueError("existing Stage-B2C seed metadata checksum mismatch")
    metadata = _read_json(metadata_path)
    structure = metadata["structure"]
    xyz_path = Path(structure["xyz"]["path"])
    cell_path = Path(structure["cell"]["path"])
    if sha256_file(xyz_path) != structure["xyz"]["sha256"]:
        raise ValueError("existing Stage-B2C coordinate checksum mismatch")
    if sha256_file(cell_path) != structure["cell"]["sha256"]:
        raise ValueError("existing Stage-B2C cell checksum mismatch")
    spec_input = analysis.get("inputs", {}).get("spec", {})
    spec_path = Path(spec_input.get("path", ""))
    if not spec_path.is_file() or sha256_file(spec_path) != spec_input.get("sha256"):
        raise ValueError("existing Stage-B2C spec checksum mismatch")
    spec = _read_json(spec_path)
    elements, positions, _comment = read_xyz(xyz_path)
    real_indices = [index for index, element in enumerate(elements) if element.upper() != "GH"]
    elements = [elements[index] for index in real_indices]
    positions = np.asarray(positions[real_indices], dtype=float)
    cell = read_cp2k_cell(cell_path)
    topology, donors = _topology_and_donors(elements, positions, cell, spec, systems, methods)
    seed_point = metadata["cavity_basis_site"]["cartesian_angstrom"]
    centroid = analysis.get("spin_density", {}).get("centroid_angstrom")
    problems: list[str] = []
    if centroid is None:
        problems.append("positive spin-density centroid is undefined")
        centroid_environment = None
    else:
        centroid_environment = _environment(centroid, donors, positions, cell, settings)
    donor_counts = Counter(int(record["molecule_id"]) for record in donors)
    spin_on_nh_bearing = 0.0
    spin_on_nh_free_amine = 0.0
    amine = spec.get("amine")
    for molecule in analysis.get("molecules", []):
        if molecule.get("component") != amine:
            continue
        fraction = float(molecule.get("positive_spin_fraction", 0.0))
        if donor_counts[int(molecule["molecule_id"])]:
            spin_on_nh_bearing += fraction
        else:
            spin_on_nh_free_amine += fraction
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "kind": "stage_b3_existing_b2c_nh_reanalysis",
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "problems": problems,
            "system_id": analysis["system_id"],
            "amine": amine,
            "replica": int(analysis["replica"]),
            "seed_role": analysis["seed_role"],
            "nh_donor_count_total": len(donors),
            "nh_donor_count_by_component": dict(
                sorted(Counter(str(record.get("component")) for record in donors).items())
            ),
            "nh_environment_at_seed": _environment(
                seed_point, donors, positions, cell, settings
            ),
            "nh_environment_at_spin_centroid": centroid_environment,
            "positive_spin_fraction_on_nh_bearing_amine_molecules": spin_on_nh_bearing,
            "positive_spin_fraction_on_nh_free_amine_molecules": spin_on_nh_free_amine,
            "vertical_attachment_proxy_ev": analysis.get("energies", {}).get(
                "vertical_attachment_proxy_ev"
            ),
            "source_analysis": {
                "path": str(Path(args.analysis).resolve()),
                "sha256": sha256_file(args.analysis),
            },
            "interpretation": settings["interpretation"],
        },
    )
    return 0


def run_reanalysis_summary(args: argparse.Namespace) -> int:
    records = [_read_json(path) for path in args.records]
    problems: list[str] = []
    keys = [
        (str(record.get("system_id")), int(record.get("replica", -1)), str(record.get("seed_role")))
        for record in records
    ]
    if len(records) != int(args.expected_count):
        problems.append(f"expected {args.expected_count} records, observed {len(records)}")
    if len(keys) != len(set(keys)):
        problems.append("duplicate system/replica/seed-role reanalysis records")
    if not records or not all(record.get("ready") for record in records):
        problems.append("one or more N-H reanalysis records are not ready")
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[str(record["system_id"])].append(record)
    systems: list[dict[str, Any]] = []
    for system_id, group in sorted(grouped.items()):
        systems.append(
            {
                "system_id": system_id,
                "amine": group[0].get("amine"),
                "record_count": len(group),
                "nh_donor_count_total": sorted(
                    {int(record["nh_donor_count_total"]) for record in group}
                ),
                "mean_positive_spin_fraction_on_nh_bearing_amine_molecules": (
                    sum(
                        float(
                            record[
                                "positive_spin_fraction_on_nh_bearing_amine_molecules"
                            ]
                        )
                        for record in group
                    )
                    / len(group)
                ),
            }
        )
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "kind": "stage_b3_existing_b2c_nh_reanalysis_summary",
            "campaign": args.campaign,
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "problems": problems,
            "record_count": len(records),
            "systems": systems,
            "records": records,
            "conclusion_scope": (
                "post hoc association diagnostics on existing fixed-nuclei PBE cubes; "
                "not a new electronic-structure calculation or hydrogen-bond proof"
            ),
        },
    )
    return 0


def run_screen(args: argparse.Namespace) -> int:
    _validate_checksum_gate(args.candidate_summary, args.candidate_gate)
    manifest = _read_json(args.candidate_manifest)
    spec = _read_json(args.spec)
    methods = _read_json(args.methods)
    systems = _read_json(args.systems)
    settings = methods["stage_b3_amine_series"]
    if str(manifest.get("system_id")) != str(spec.get("system_id")) or int(
        manifest.get("replica", -1)
    ) != int(spec.get("replica", -2)):
        raise ValueError("Stage-B3 source candidate manifest and spec differ")
    elements, positions, cell_path, cell, structure = _source_structure(manifest)
    topology, donors = _topology_and_donors(elements, positions, cell, spec, systems, methods)
    sites = list(manifest.get("void_search", {}).get("sites", []))
    selected = select_nh_contrast_sites(
        sites,
        donors,
        positions,
        cell,
        shell_radii_angstrom=settings["nh_shell_radii_angstrom"],
        selection_radius_angstrom=float(settings["nh_selection_radius_angstrom"]),
        inward_angle_minimum_degrees=float(settings["inward_angle_minimum_degrees"]),
        distance_decay_angstrom=float(settings["distance_decay_angstrom"]),
        minimum_separation_angstrom=float(settings["minimum_seed_separation_angstrom"]),
    )
    centers = molecular_heavy_atom_centers(elements, positions, cell, topology)
    counts = _component_counts_from_spec(spec)
    records: list[dict[str, Any]] = []
    for record in selected:
        records.append(
            {
                **record,
                "local_composition_at_seed": local_composition(
                    record["site"]["cartesian_angstrom"],
                    centers,
                    cell,
                    settings["composition_shell_radii_angstrom"],
                    counts,
                    amine_component=spec.get("amine"),
                ),
            }
        )
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "kind": "stage_b3_nh_environment_screen",
            "ready": len(records) == 2,
            "scientific_status": SCIENTIFIC_STATUS,
            "campaign": args.campaign,
            "source_campaign": args.source_campaign,
            "system_id": spec["system_id"],
            "amine": spec.get("amine"),
            "replica": int(spec["replica"]),
            "component_counts": counts,
            "nh_donor_count_total": len(donors),
            "nh_donor_count_by_component": dict(
                sorted(Counter(str(record.get("component")) for record in donors).items())
            ),
            "selected_sites": records,
            "source": {
                "candidate_manifest": {
                    "path": str(Path(args.candidate_manifest).resolve()),
                    "sha256": sha256_file(args.candidate_manifest),
                },
                "spec": {"path": str(Path(args.spec).resolve()), "sha256": sha256_file(args.spec)},
                "xyz": structure["xyz"],
                "cell": {"path": str(cell_path), "sha256": sha256_file(cell_path)},
            },
            "interpretation": settings["interpretation"],
        },
    )
    return 0


def run_environment_summary(args: argparse.Namespace) -> int:
    records = [_read_json(path) for path in args.records]
    expected = {str(value) for value in args.expected_systems}
    observed = {str(record.get("system_id")) for record in records}
    problems: list[str] = []
    if observed != expected:
        problems.append(f"environment systems differ from expected: {sorted(observed)}")
    if not records or not all(record.get("ready") for record in records):
        problems.append("one or more amine environment screens are not ready")
    table = []
    for record in sorted(records, key=lambda value: str(value["system_id"])):
        sites = {value["seed_role"]: value for value in record["selected_sites"]}
        table.append(
            {
                "system_id": record["system_id"],
                "amine": record.get("amine"),
                "nh_donor_count_total": record["nh_donor_count_total"],
                "nh_degenerate_reference": bool(sites["nh_facing"]["nh_degenerate_reference"]),
                "selection_inward_nh_count": {
                    role: sites[role]["selection_inward_nh_count"]
                    for role in ("nh_facing", "nh_control")
                },
                "selection_ranking_score": {
                    role: sites[role]["selection_ranking_score"]
                    for role in ("nh_facing", "nh_control")
                },
            }
        )
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "kind": "stage_b3_amine_environment_summary",
            "campaign": args.campaign,
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "problems": problems,
            "systems": table,
            "records": records,
        },
    )
    return 0


def run_prepare(args: argparse.Namespace) -> int:
    screen = _read_json(args.screen)
    manifest = _read_json(args.candidate_manifest)
    if not screen.get("ready"):
        raise ValueError("Stage-B3 N-H environment screen is not ready")
    elements, positions, cell_path, _cell, _structure = _source_structure(manifest)
    output_dir = Path(args.output_dir).resolve()
    records: list[dict[str, Any]] = []
    for selected in screen["selected_sites"]:
        role = str(selected["seed_role"])
        seed_dir = output_dir / role
        xyz = seed_dir / "coordinates.xyz"
        cell_output = seed_dir / "cell.inc"
        metadata = seed_dir / "metadata.json"
        center = np.asarray(selected["site"]["cartesian_angstrom"], dtype=float)
        write_xyz(
            xyz,
            [*elements, "Gh"],
            np.vstack([positions, center]),
            (
                f"stage_b3_seed={role} source={screen['system_id']}_r{screen['replica']} "
                "Li_absent PFAS_absent fixed_nuclei Gh=basis_seed_not_a_restraint"
            ),
        )
        cell_output.parent.mkdir(parents=True, exist_ok=True)
        cell_output.write_text(cell_path.read_text(encoding="utf-8"), encoding="utf-8")
        record = {
            "schema_version": 1,
            "ready": True,
            "scientific_status": SCIENTIFIC_STATUS,
            "system_id": screen["system_id"],
            "amine": screen.get("amine"),
            "replica": int(screen["replica"]),
            **selected,
            "structure": {
                "solvent_atom_count": len(elements),
                "ghost_basis_center_count": 1,
                "xyz": {"path": str(xyz), "sha256": sha256_file(xyz)},
                "cell": {"path": str(cell_output), "sha256": sha256_file(cell_output)},
            },
            "source_screen": {
                "path": str(Path(args.screen).resolve()),
                "sha256": sha256_file(args.screen),
            },
        }
        _write_json(metadata, record)
        record["metadata"] = {"path": str(metadata), "sha256": sha256_file(metadata)}
        records.append(record)
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "kind": "stage_b3_electronic_seed_manifest",
            "ready": len(records) == 2,
            "scientific_status": SCIENTIFIC_STATUS,
            "system_id": screen["system_id"],
            "amine": screen.get("amine"),
            "replica": int(screen["replica"]),
            "seed_roles": [record["seed_role"] for record in records],
            "seeds": records,
        },
    )
    return 0


def _ready_seed(manifest: dict[str, Any], role: str) -> dict[str, Any]:
    records = [
        record for record in manifest.get("seeds", []) if record.get("seed_role") == role
    ]
    if not manifest.get("ready") or len(records) != 1 or not records[0].get("ready"):
        raise ValueError(f"Stage-B3 seed role {role!r} is not uniquely ready")
    return records[0]


def _states(settings: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {str(record["id"]): record for record in settings["states"]}


def run_render(args: argparse.Namespace) -> int:
    manifest = _read_json(args.manifest)
    settings = _read_json(args.methods)["stage_b3_amine_series"]
    seed = _ready_seed(manifest, args.seed_role)
    structure = seed["structure"]
    xyz = Path(structure["xyz"]["path"])
    cell = Path(structure["cell"]["path"])
    if sha256_file(xyz) != structure["xyz"]["sha256"]:
        raise ValueError("Stage-B3 seed coordinate checksum mismatch")
    if sha256_file(cell) != structure["cell"]["sha256"]:
        raise ValueError("Stage-B3 seed cell checksum mismatch")
    for state_id, output in (("neutral", args.neutral_output), ("anion", args.anion_output)):
        render_stage_b3_amine_series_cp2k(
            args.template,
            output,
            project=f"{args.project_prefix}_{state_id}_stage_b3",
            coordinates_path=xyz,
            cell_path=cell,
            method=settings,
            state=_states(settings)[state_id],
        )
    return 0


def _input_problems(text: str, state: dict[str, Any]) -> list[str]:
    problems: list[str] = []
    charge = int(state["charge"])
    multiplicity = int(state["multiplicity"])
    expected_uks = "TRUE" if state["uks"] else "FALSE"
    if not re.search(rf"^\s*CHARGE\s+{charge}\s*$", text, flags=re.MULTILINE):
        problems.append(f"CP2K input is not explicit charge {charge}")
    if not re.search(rf"^\s*MULTIPLICITY\s+{multiplicity}\s*$", text, flags=re.MULTILINE):
        problems.append(f"CP2K input is not explicit multiplicity {multiplicity}")
    if not re.search(rf"^\s*UKS\s+{expected_uks}\s*$", text, flags=re.MULTILINE):
        problems.append(f"CP2K input does not set UKS {expected_uks}")
    upper = text.upper()
    if "RUN_TYPE ENERGY" not in upper:
        problems.append("CP2K input is not a fixed-nuclei ENERGY calculation")
    if "&CDFT" in upper or "&CONSTRAINT" in upper or "&KIND LI" in upper:
        problems.append("CP2K input contains forbidden Li/localization machinery")
    return problems


def run_analyze(args: argparse.Namespace) -> int:
    manifest = _read_json(args.manifest)
    seed = _ready_seed(manifest, args.seed_role)
    spec = _read_json(args.spec)
    methods = _read_json(args.methods)
    systems = _read_json(args.systems)
    settings = methods["stage_b3_amine_series"]
    states = _states(settings)
    structure = seed["structure"]
    xyz = Path(structure["xyz"]["path"])
    cell_path = Path(structure["cell"]["path"])
    elements, positions, _comment = read_xyz(xyz)
    positions = np.asarray(positions, dtype=float)
    cell = read_cp2k_cell(cell_path)
    topology, donors = _topology_and_donors(
        elements, positions, cell, spec, systems, methods
    )
    problems: list[str] = []
    neutral_input = Path(args.neutral_input)
    anion_input = Path(args.anion_input)
    neutral_output = Path(args.neutral_output)
    anion_output = Path(args.anion_output)
    problems.extend(_input_problems(neutral_input.read_text(encoding="utf-8"), states["neutral"]))
    problems.extend(_input_problems(anion_input.read_text(encoding="utf-8"), states["anion"]))
    neutral = parse_cp2k_text(neutral_output.read_text(encoding="utf-8", errors="replace"))
    anion = parse_cp2k_text(anion_output.read_text(encoding="utf-8", errors="replace"))
    problems.extend(f"neutral: {problem}" for problem in neutral.problems)
    problems.extend(f"anion: {problem}" for problem in anion.problems)
    spin_path = Path(args.spin_cube)
    electron_path = Path(args.electron_cube)
    spin_cube = read_cube(spin_path)
    electron_cube = read_cube(electron_path)
    tolerance = float(methods["stage_b_localization"]["cube_grid_tolerance_angstrom"])
    problems.extend(cube_compatibility_problems(spin_cube, electron_cube, tolerance_angstrom=tolerance))
    problems.extend(
        candidate_cube_problems(
            spin_cube, elements, np.asarray(positions, dtype=float), cell, tolerance_angstrom=tolerance
        )
    )
    seed_center = np.asarray(seed["site"]["cartesian_angstrom"], dtype=float)
    localization_settings = {
        **methods["stage_b_localization"],
        "cavity_probe_radius_angstrom": float(settings["cavity_probe_radius_angstrom"]),
    }
    localization = state_localization_record(
        spin_cube,
        elements,
        np.asarray(positions, dtype=float),
        cell,
        topology,
        seed_center,
        localization_settings,
    )
    signed = float(localization["spin_density"]["signed_integral"])
    if not float(localization_settings["signed_spin_integral_min"]) <= signed <= float(
        localization_settings["signed_spin_integral_max"]
    ):
        problems.append(f"signed spin integral {signed:.6g} is outside the accepted range")
    centroid = localization["spin_density"].get("centroid_angstrom")
    centroid_environment = None
    if centroid is None:
        problems.append("positive spin-density centroid is undefined")
    else:
        centroid_environment = _environment(
            centroid, donors, positions, cell, settings
        )
    component_spin: dict[str, float] = defaultdict(float)
    for molecule in localization.get("molecules", []):
        component_spin[str(molecule.get("component"))] += float(
            molecule["positive_spin_fraction"]
        )
    attachment_hartree = None
    attachment_ev = None
    if neutral.energy_hartree is not None and anion.energy_hartree is not None:
        attachment_hartree = neutral.energy_hartree - anion.energy_hartree
        attachment_ev = attachment_hartree * HARTREE_TO_EV
    else:
        problems.append("vertical attachment-energy proxy could not be evaluated")
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "kind": "stage_b3_electronic_pair",
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "problems": problems,
            "system_id": manifest["system_id"],
            "amine": manifest.get("amine"),
            "replica": int(manifest["replica"]),
            "seed_role": args.seed_role,
            "nh_degenerate_reference": bool(seed["nh_degenerate_reference"]),
            "energies": {
                "neutral_hartree": neutral.energy_hartree,
                "anion_hartree": anion.energy_hartree,
                "vertical_attachment_proxy_hartree": attachment_hartree,
                "vertical_attachment_proxy_ev": attachment_ev,
                "definition": "E(neutral)-E(anion) at identical nuclei and ghost basis seed",
                "finite_size_corrected": False,
            },
            "nh_environment_at_seed": seed["nh_environment"],
            "nh_environment_at_spin_centroid": centroid_environment,
            "local_composition_at_seed": seed["local_composition_at_seed"],
            "positive_spin_fraction_by_component": dict(sorted(component_spin.items())),
            "inputs": {
                "neutral_input": {"path": str(neutral_input.resolve()), "sha256": sha256_file(neutral_input)},
                "neutral_output": {"path": str(neutral_output.resolve()), "sha256": sha256_file(neutral_output)},
                "anion_input": {"path": str(anion_input.resolve()), "sha256": sha256_file(anion_input)},
                "anion_output": {"path": str(anion_output.resolve()), "sha256": sha256_file(anion_output)},
                "spin_cube": {"path": str(spin_path.resolve()), "sha256": sha256_file(spin_path)},
                "electron_cube": {"path": str(electron_path.resolve()), "sha256": sha256_file(electron_path)},
                "seed_metadata": seed["metadata"],
                "spec": {"path": str(Path(args.spec).resolve()), "sha256": sha256_file(args.spec)},
            },
            **localization,
        },
    )
    return 0


def run_pilot_summary(args: argparse.Namespace) -> int:
    records = [_read_json(path) for path in args.records]
    expected = {
        (system, int(replica), role)
        for system in args.expected_systems
        for replica in args.expected_replicas
        for role in ("nh_facing", "nh_control")
    }
    observed = {
        (str(record.get("system_id")), int(record.get("replica", -1)), str(record.get("seed_role")))
        for record in records
    }
    problems: list[str] = []
    if observed != expected:
        problems.append(f"pilot keys differ from expected: {len(observed)} versus {len(expected)}")
    if len(records) != len(observed):
        problems.append("duplicate Stage-B3 pilot records")
    if not records or not all(record.get("ready") for record in records):
        problems.append("one or more Stage-B3 electronic records are not ready")
    grouped: dict[tuple[str, int], dict[str, dict[str, Any]]] = defaultdict(dict)
    for record in records:
        grouped[(str(record["system_id"]), int(record["replica"]))][
            str(record["seed_role"])
        ] = record
    comparisons: list[dict[str, Any]] = []
    tolerance = float(args.energy_tolerance_ev)
    for (system_id, replica), roles in sorted(grouped.items()):
        if set(roles) != {"nh_facing", "nh_control"}:
            continue
        facing = roles["nh_facing"]
        control = roles["nh_control"]
        facing_energy = facing["energies"]["vertical_attachment_proxy_ev"]
        control_energy = control["energies"]["vertical_attachment_proxy_ev"]
        delta = (
            float(facing_energy) - float(control_energy)
            if facing_energy is not None and control_energy is not None
            else None
        )
        degenerate = bool(facing["nh_degenerate_reference"] and control["nh_degenerate_reference"])
        if degenerate:
            result = "nh_free_or_geometrically_degenerate_reference"
        elif delta is None:
            result = "unavailable"
        elif abs(delta) <= tolerance:
            result = "indistinguishable_at_pilot_tolerance"
        elif delta > 0:
            result = "nh_facing_has_larger_attachment_proxy"
        else:
            result = "nh_control_has_larger_attachment_proxy"
        comparisons.append(
            {
                "system_id": system_id,
                "amine": facing.get("amine"),
                "replica": replica,
                "nh_degenerate_reference": degenerate,
                "vertical_attachment_proxy_ev": {
                    "nh_facing": facing_energy,
                    "nh_control": control_energy,
                    "facing_minus_control": delta,
                },
                "pilot_result": result,
                "energy_tolerance_ev": tolerance,
            }
        )
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "kind": "stage_b3_amine_series_electronic_pilot",
            "campaign": args.campaign,
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "problems": problems,
            "record_count": len(records),
            "comparisons": comparisons,
            "records": records,
            "method_limitations": [
                "fixed nuclei sample pre-existing solvent fluctuations only",
                "PBE may over-delocalize the excess electron",
                "charged periodic cells use a compensating background and are not finite-size corrected",
                "N-H-facing is a geometric association label, not proof of a hydrogen bond",
                "Li and PFAS are absent; this is not Li ionization thermodynamics or reaction kinetics",
            ],
        },
    )
    return 0


def run_combine(args: argparse.Namespace) -> int:
    for summary, gate in (
        (args.reanalysis_summary, args.reanalysis_gate),
        (args.environment_summary, args.environment_gate),
        (args.pilot_summary, args.pilot_gate),
    ):
        _validate_checksum_gate(summary, gate)
    summaries = {
        "existing_b2c_nh_reanalysis": _read_json(args.reanalysis_summary),
        "amine_environment_panel": _read_json(args.environment_summary),
        "electronic_pilot": _read_json(args.pilot_summary),
    }
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "kind": "stage_b3_amine_series",
            "campaign": args.campaign,
            "ready": all(value.get("ready") for value in summaries.values()),
            "scientific_status": SCIENTIFIC_STATUS,
            "summaries": summaries,
            "next_gate": (
                "Only after this screen: choose the smallest representative subset for a "
                "range-separated/hybrid and finite-size sensitivity check before paper claims."
            ),
        },
    )
    return 0


def run_gate(args: argparse.Namespace) -> int:
    summary = Path(args.summary).resolve()
    if not _read_json(summary).get("ready"):
        raise ValueError(f"Stage-B3 summary is not ready: {summary}")
    digest = hashlib.sha256(summary.read_bytes()).hexdigest()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(f"sha256 {digest}  {summary.name}\n", encoding="utf-8")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    reanalyze = subparsers.add_parser("reanalyze")
    reanalyze.add_argument("--analysis", required=True)
    reanalyze.add_argument("--methods", required=True)
    reanalyze.add_argument("--systems", required=True)
    reanalyze.add_argument("--output", required=True)
    reanalyze.set_defaults(func=run_reanalyze)

    reanalysis_summary = subparsers.add_parser("reanalysis-summary")
    reanalysis_summary.add_argument("--campaign", required=True)
    reanalysis_summary.add_argument("--expected-count", type=int, required=True)
    reanalysis_summary.add_argument("--output", required=True)
    reanalysis_summary.add_argument("--records", nargs="+", required=True)
    reanalysis_summary.set_defaults(func=run_reanalysis_summary)

    screen = subparsers.add_parser("screen")
    screen.add_argument("--campaign", required=True)
    screen.add_argument("--source-campaign", required=True)
    screen.add_argument("--candidate-manifest", required=True)
    screen.add_argument("--candidate-summary", required=True)
    screen.add_argument("--candidate-gate", required=True)
    screen.add_argument("--spec", required=True)
    screen.add_argument("--methods", required=True)
    screen.add_argument("--systems", required=True)
    screen.add_argument("--output", required=True)
    screen.set_defaults(func=run_screen)

    environment_summary = subparsers.add_parser("environment-summary")
    environment_summary.add_argument("--campaign", required=True)
    environment_summary.add_argument("--expected-systems", nargs="+", required=True)
    environment_summary.add_argument("--output", required=True)
    environment_summary.add_argument("--records", nargs="+", required=True)
    environment_summary.set_defaults(func=run_environment_summary)

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--screen", required=True)
    prepare.add_argument("--candidate-manifest", required=True)
    prepare.add_argument("--output-dir", required=True)
    prepare.add_argument("--output", required=True)
    prepare.set_defaults(func=run_prepare)

    render = subparsers.add_parser("render")
    render.add_argument("--manifest", required=True)
    render.add_argument("--methods", required=True)
    render.add_argument("--template", required=True)
    render.add_argument("--seed-role", required=True)
    render.add_argument("--project-prefix", required=True)
    render.add_argument("--neutral-output", required=True)
    render.add_argument("--anion-output", required=True)
    render.set_defaults(func=run_render)

    analyze = subparsers.add_parser("analyze")
    analyze.add_argument("--manifest", required=True)
    analyze.add_argument("--seed-role", required=True)
    analyze.add_argument("--spec", required=True)
    analyze.add_argument("--methods", required=True)
    analyze.add_argument("--systems", required=True)
    analyze.add_argument("--neutral-input", required=True)
    analyze.add_argument("--neutral-output", required=True)
    analyze.add_argument("--anion-input", required=True)
    analyze.add_argument("--anion-output", required=True)
    analyze.add_argument("--spin-cube", required=True)
    analyze.add_argument("--electron-cube", required=True)
    analyze.add_argument("--output", required=True)
    analyze.set_defaults(func=run_analyze)

    pilot_summary = subparsers.add_parser("pilot-summary")
    pilot_summary.add_argument("--campaign", required=True)
    pilot_summary.add_argument("--expected-systems", nargs="+", required=True)
    pilot_summary.add_argument("--expected-replicas", nargs="+", type=int, required=True)
    pilot_summary.add_argument("--energy-tolerance-ev", type=float, required=True)
    pilot_summary.add_argument("--output", required=True)
    pilot_summary.add_argument("--records", nargs="+", required=True)
    pilot_summary.set_defaults(func=run_pilot_summary)

    combine = subparsers.add_parser("combine")
    combine.add_argument("--campaign", required=True)
    combine.add_argument("--reanalysis-summary", required=True)
    combine.add_argument("--reanalysis-gate", required=True)
    combine.add_argument("--environment-summary", required=True)
    combine.add_argument("--environment-gate", required=True)
    combine.add_argument("--pilot-summary", required=True)
    combine.add_argument("--pilot-gate", required=True)
    combine.add_argument("--output", required=True)
    combine.set_defaults(func=run_combine)

    gate = subparsers.add_parser("gate")
    gate.add_argument("--summary", required=True)
    gate.add_argument("--output", required=True)
    gate.set_defaults(func=run_gate)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        return int(args.func(args))
    except (KeyError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
