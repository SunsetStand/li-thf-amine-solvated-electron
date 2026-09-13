#!/usr/bin/env python3
"""Prepare and analyze fixed-nuclei Stage-B2C preferential-solvation smoke pairs."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import defaultdict
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
from solvelec.parsers import HARTREE_TO_EV, parse_cp2k_text
from solvelec.preferential_solvation import (
    local_composition,
    molecular_heavy_atom_centers,
    select_preferential_sites,
)
from solvelec.provenance import sha256_file
from solvelec.rendering import render_stage_b2c_preferential_cp2k

SCIENTIFIC_STATUS = "NUMERICAL_PREFERENTIAL_SOLVATION_SMOKE_ONLY"


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


def _validate_checksum_gate(summary_path: str | Path, gate_path: str | Path) -> None:
    summary = Path(summary_path)
    gate = Path(gate_path)
    if not summary.is_file() or not gate.is_file():
        raise ValueError("Stage-B2C requires an accepted Stage-B candidate summary and gate")
    if not _read_json(summary).get("ready"):
        raise ValueError("Stage-B2C source candidate summary is not ready")
    digest = hashlib.sha256(summary.read_bytes()).hexdigest()
    if gate.read_text(encoding="utf-8").split() != ["sha256", digest, summary.name]:
        raise ValueError("Stage-B2C source candidate gate does not match its summary")


def _component_counts_from_spec(spec: dict[str, Any]) -> dict[str, int]:
    configured = spec.get("component_counts")
    if isinstance(configured, dict) and configured:
        counts = {str(name): int(value) for name, value in configured.items()}
    else:
        counts: dict[str, int] = {}
        thf_count = int(spec.get("thf_count", 0))
        if thf_count:
            counts["thf"] = thf_count
        amine = spec.get("amine")
        amine_count = int(spec.get("amine_count_initial", 0))
        if amine and amine_count:
            counts[str(amine)] = amine_count
    if not counts or any(value <= 0 for value in counts.values()):
        raise ValueError("Stage-B2C spec has no positive solvent component count")
    return counts


def _component_definitions(systems: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        "thf": systems["thf"],
        **{str(name): value for name, value in systems.get("amines", {}).items()},
    }


def _ready_seed(manifest: dict[str, Any], seed_role: str) -> dict[str, Any]:
    records = [
        record
        for record in manifest.get("seeds", [])
        if record.get("seed_role") == seed_role and record.get("ready")
    ]
    if not manifest.get("ready") or len(records) != 1:
        raise ValueError(f"Stage-B2C seed role {seed_role!r} is not uniquely ready")
    return records[0]


def _states(settings: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {str(record["id"]): record for record in settings["states"]}


def run_prepare(args: argparse.Namespace) -> int:
    _validate_checksum_gate(args.candidate_summary, args.candidate_gate)
    candidate_manifest_path = Path(args.candidate_manifest).resolve()
    spec_path = Path(args.spec).resolve()
    methods_path = Path(args.methods).resolve()
    systems_path = Path(args.systems).resolve()
    candidate_manifest = _read_json(candidate_manifest_path)
    spec = _read_json(spec_path)
    methods = _read_json(methods_path)
    systems = _read_json(systems_path)
    settings = methods["stage_b2c_preferential_smoke"]
    if settings.get("scientific_status") != SCIENTIFIC_STATUS:
        raise ValueError("Stage-B2C method is not explicitly smoke-only")
    if not candidate_manifest.get("ready"):
        raise ValueError("Stage-B2C source candidate manifest is not ready")
    if str(candidate_manifest.get("system_id")) != str(spec.get("system_id")) or int(
        candidate_manifest.get("replica", -1)
    ) != int(spec.get("replica", -2)):
        raise ValueError("Stage-B2C source candidate manifest and solvent spec differ")

    candidates = [
        record for record in candidate_manifest.get("candidates", []) if record.get("ready")
    ]
    if not candidates:
        raise ValueError("Stage-B2C source candidate manifest has no ready structure")
    structure = candidates[0]["structure"]
    source_xyz = Path(structure["xyz"]["path"])
    source_cell = Path(structure["cell"]["path"])
    if sha256_file(source_xyz) != structure["xyz"]["sha256"]:
        raise ValueError("Stage-B2C source candidate XYZ checksum mismatch")
    if sha256_file(source_cell) != structure["cell"]["sha256"]:
        raise ValueError("Stage-B2C source candidate cell checksum mismatch")
    elements, positions, _comment = read_xyz(source_xyz)
    li_indices = [index for index, element in enumerate(elements) if element.upper() == "LI"]
    ghost_indices = [index for index, element in enumerate(elements) if element.upper() == "GH"]
    if len(li_indices) != 1 or len(ghost_indices) != 1:
        raise ValueError("Stage-B2C source must contain exactly one Li and one ghost")
    solvent_indices = [
        index for index, element in enumerate(elements) if element.upper() not in {"LI", "GH"}
    ]
    solvent_elements = [elements[index] for index in solvent_indices]
    solvent_positions = positions[solvent_indices]
    cell = read_cp2k_cell(source_cell)
    component_counts = _component_counts_from_spec(spec)
    topology = infer_molecular_topology(
        solvent_elements,
        solvent_positions,
        cell,
        component_counts,
        _component_definitions(systems),
        bond_scale=float(methods["stage_b_localization"]["covalent_bond_scale"]),
    )
    if topology.problems:
        raise ValueError("Stage-B2C molecular topology failed: " + "; ".join(topology.problems))
    molecule_centers = molecular_heavy_atom_centers(
        solvent_elements, solvent_positions, cell, topology
    )
    sites = list(candidate_manifest.get("void_search", {}).get("sites", []))
    selected = select_preferential_sites(
        sites,
        molecule_centers,
        cell,
        component_counts,
        amine_component=spec.get("amine"),
        shell_radii_angstrom=settings["local_shell_radii_angstrom"],
        selection_radius_angstrom=float(settings["selection_radius_angstrom"]),
        minimum_separation_angstrom=float(settings["minimum_seed_separation_angstrom"]),
    )
    output_dir = Path(args.output_dir).resolve()
    records: list[dict[str, Any]] = []
    for selection in selected:
        role = str(selection["seed_role"])
        seed_dir = output_dir / role
        coordinates = seed_dir / "coordinates.xyz"
        cell_output = seed_dir / "cell.inc"
        metadata = seed_dir / "metadata.json"
        seed_position = np.asarray(selection["site"]["cartesian_angstrom"], dtype=float)
        write_xyz(
            coordinates,
            [*solvent_elements, "Gh"],
            np.vstack([solvent_positions, seed_position]),
            (
                f"stage_b2c_seed={role} source={candidate_manifest['snapshot_id']} "
                "Li_absent PFAS_absent fixed_nuclei Gh=basis_seed_not_a_restraint"
            ),
        )
        cell_output.parent.mkdir(parents=True, exist_ok=True)
        cell_output.write_text(source_cell.read_text(encoding="utf-8"), encoding="utf-8")
        record = {
            "schema_version": 1,
            "ready": True,
            "scientific_status": SCIENTIFIC_STATUS,
            "campaign": args.campaign,
            "source_campaign": args.source_campaign,
            "system_id": spec["system_id"],
            "amine": spec.get("amine"),
            "replica": int(spec["replica"]),
            "snapshot_id": candidate_manifest["snapshot_id"],
            "seed_role": role,
            "void_rank": int(selection["site"]["rank"]),
            "cavity_basis_site": selection["site"],
            "local_composition_at_seed": selection["local_composition"],
            "selection_amine_mole_fraction": selection["selection_amine_mole_fraction"],
            "composition_degenerate_reference": selection["composition_degenerate_reference"],
            "electronic_pair": {
                "states": settings["states"],
                "li_atom_present": False,
                "pfas_present": False,
                "localization_constraint_present": False,
                "fixed_nuclei": True,
                "same_geometry_and_basis_seed": True,
            },
            "structure": {
                "solvent_atom_count": len(solvent_elements),
                "atom_count_including_ghost": len(solvent_elements) + 1,
                "ghost_basis_center_count": 1,
                "component_counts": component_counts,
                "xyz": {"path": str(coordinates), "sha256": sha256_file(coordinates)},
                "cell": {"path": str(cell_output), "sha256": sha256_file(cell_output)},
            },
            "source": {
                "candidate_manifest": {
                    "path": str(candidate_manifest_path),
                    "sha256": sha256_file(candidate_manifest_path),
                },
                "candidate_summary": {
                    "path": str(Path(args.candidate_summary).resolve()),
                    "sha256": sha256_file(args.candidate_summary),
                },
                "candidate_gate": {
                    "path": str(Path(args.candidate_gate).resolve()),
                    "sha256": sha256_file(args.candidate_gate),
                },
                "spec": {"path": str(spec_path), "sha256": sha256_file(spec_path)},
                "methods": {"path": str(methods_path), "sha256": sha256_file(methods_path)},
                "systems": {"path": str(systems_path), "sha256": sha256_file(systems_path)},
            },
        }
        _write_json(metadata, record)
        record["metadata"] = {"path": str(metadata), "sha256": sha256_file(metadata)}
        records.append(record)

    output = {
        "schema_version": 1,
        "ready": len(records) == 2,
        "scientific_status": SCIENTIFIC_STATUS,
        "campaign": args.campaign,
        "source_campaign": args.source_campaign,
        "system_id": spec["system_id"],
        "amine": spec.get("amine"),
        "replica": int(spec["replica"]),
        "snapshot_id": candidate_manifest["snapshot_id"],
        "seed_roles": [record["seed_role"] for record in records],
        "component_counts": component_counts,
        "molecule_centers": list(molecule_centers),
        "seeds": records,
    }
    _write_json(args.output, output)
    return 0


def run_render(args: argparse.Namespace) -> int:
    manifest = _read_json(args.manifest)
    settings = _read_json(args.methods)["stage_b2c_preferential_smoke"]
    seed = _ready_seed(manifest, args.seed_role)
    structure = seed["structure"]
    coordinates = Path(structure["xyz"]["path"])
    cell = Path(structure["cell"]["path"])
    if sha256_file(coordinates) != structure["xyz"]["sha256"]:
        raise ValueError("Stage-B2C coordinate checksum mismatch")
    if sha256_file(cell) != structure["cell"]["sha256"]:
        raise ValueError("Stage-B2C cell checksum mismatch")
    elements, _positions, _comment = read_xyz(coordinates)
    if any(element.upper() == "LI" for element in elements):
        raise ValueError("Stage-B2C coordinates unexpectedly contain Li")
    states = _states(settings)
    for state_id, output in (("neutral", args.neutral_output), ("anion", args.anion_output)):
        render_stage_b2c_preferential_cp2k(
            args.template,
            output,
            project=f"{args.project_prefix}_{state_id}_stage_b2c_preferential_smoke",
            coordinates_path=coordinates,
            cell_path=cell,
            method=settings,
            state=states[state_id],
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
    if "&CDFT" in upper or "&CONSTRAINT" in upper:
        problems.append("CP2K input contains a forbidden localization constraint")
    if "&KIND LI" in upper:
        problems.append("CP2K input contains a forbidden Li kind")
    return problems


def run_analyze(args: argparse.Namespace) -> int:
    manifest = _read_json(args.manifest)
    seed = _ready_seed(manifest, args.seed_role)
    spec = _read_json(args.spec)
    methods = _read_json(args.methods)
    systems = _read_json(args.systems)
    settings = methods["stage_b2c_preferential_smoke"]
    states = _states(settings)
    if str(seed["system_id"]) != args.system or int(seed["replica"]) != int(args.replica):
        raise ValueError("Stage-B2C seed does not match requested system/replica")
    if str(spec.get("system_id")) != args.system or int(spec.get("replica", -1)) != int(
        args.replica
    ):
        raise ValueError("Stage-B2C spec does not match requested system/replica")
    structure = seed["structure"]
    coordinates = Path(structure["xyz"]["path"])
    cell_path = Path(structure["cell"]["path"])
    for path, expected in (
        (coordinates, structure["xyz"]["sha256"]),
        (cell_path, structure["cell"]["sha256"]),
    ):
        if sha256_file(path) != expected:
            raise ValueError(f"Stage-B2C immutable structure checksum mismatch: {path}")
    elements, positions, _comment = read_xyz(coordinates)
    if any(element.upper() == "LI" for element in elements):
        raise ValueError("Stage-B2C analysis refuses a structure containing Li")
    if sum(element.upper() == "GH" for element in elements) != 1:
        raise ValueError("Stage-B2C structure must contain exactly one ghost center")
    cell = read_cp2k_cell(cell_path)
    localization_settings = {
        **methods["stage_b_localization"],
        "cavity_probe_radius_angstrom": float(settings["cavity_probe_radius_angstrom"]),
    }
    topology = infer_molecular_topology(
        elements,
        positions,
        cell,
        _component_counts_from_spec(spec),
        _component_definitions(systems),
        bond_scale=float(localization_settings["covalent_bond_scale"]),
    )
    problems = [*topology.problems]
    neutral_input = Path(args.neutral_input)
    anion_input = Path(args.anion_input)
    neutral_output = Path(args.neutral_output)
    anion_output = Path(args.anion_output)
    problems.extend(_input_problems(neutral_input.read_text(encoding="utf-8"), states["neutral"]))
    problems.extend(_input_problems(anion_input.read_text(encoding="utf-8"), states["anion"]))
    neutral_result = parse_cp2k_text(neutral_output.read_text(encoding="utf-8", errors="replace"))
    anion_result = parse_cp2k_text(anion_output.read_text(encoding="utf-8", errors="replace"))
    problems.extend(f"neutral: {problem}" for problem in neutral_result.problems)
    problems.extend(f"anion: {problem}" for problem in anion_result.problems)
    spin_path = Path(args.spin_cube)
    electron_path = Path(args.electron_cube)
    spin_cube = read_cube(spin_path)
    electron_cube = read_cube(electron_path)
    tolerance = float(localization_settings["cube_grid_tolerance_angstrom"])
    problems.extend(
        cube_compatibility_problems(spin_cube, electron_cube, tolerance_angstrom=tolerance)
    )
    problems.extend(
        candidate_cube_problems(spin_cube, elements, positions, cell, tolerance_angstrom=tolerance)
    )
    seed_center = np.asarray(seed["cavity_basis_site"]["cartesian_angstrom"], dtype=float)
    localization = state_localization_record(
        spin_cube, elements, positions, cell, topology, seed_center, localization_settings
    )
    signed = float(localization["spin_density"]["signed_integral"])
    lower = float(localization_settings["signed_spin_integral_min"])
    upper = float(localization_settings["signed_spin_integral_max"])
    if not lower <= signed <= upper:
        problems.append(f"signed spin integral {signed:.6g} is outside [{lower:.6g}, {upper:.6g}]")
    centroid = localization["spin_density"].get("centroid_angstrom")
    final_composition = None
    if centroid is None:
        problems.append("positive spin-density centroid is undefined")
    else:
        centers = molecular_heavy_atom_centers(elements, positions, cell, topology)
        final_composition = local_composition(
            centroid,
            centers,
            cell,
            settings["local_shell_radii_angstrom"],
            _component_counts_from_spec(spec),
            amine_component=spec.get("amine"),
        )
    component_spin: dict[str, float] = defaultdict(float)
    for molecule in localization.get("molecules", []):
        component_spin[str(molecule.get("component"))] += float(
            molecule["positive_spin_fraction"]
        )
    energy_hartree = None
    energy_ev = None
    if neutral_result.energy_hartree is not None and anion_result.energy_hartree is not None:
        energy_hartree = neutral_result.energy_hartree - anion_result.energy_hartree
        energy_ev = energy_hartree * HARTREE_TO_EV
    else:
        problems.append("vertical attachment-energy proxy could not be evaluated")
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "campaign": args.campaign,
            "source_campaign": seed["source_campaign"],
            "kind": "stage_b2c_preferential_pair",
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "system_id": args.system,
            "amine": spec.get("amine"),
            "replica": int(args.replica),
            "seed_role": args.seed_role,
            "void_rank": int(seed["void_rank"]),
            "problems": problems,
            "checks": {
                "li_absent": True,
                "pfas_absent": True,
                "fixed_nuclei_pair": not any("fixed-nuclei" in problem for problem in problems),
                "constraint_release": not any("constraint" in problem for problem in problems),
                "molecular_topology": not topology.problems,
                "signed_spin_integral": not any(
                    problem.startswith("signed spin integral") for problem in problems
                ),
            },
            "energies": {
                "neutral_hartree": neutral_result.energy_hartree,
                "anion_hartree": anion_result.energy_hartree,
                "vertical_attachment_proxy_hartree": energy_hartree,
                "vertical_attachment_proxy_ev": energy_ev,
                "definition": "E(neutral)-E(anion) at identical nuclei and ghost basis seed",
                "finite_size_corrected": False,
            },
            "local_composition_at_seed": seed["local_composition_at_seed"],
            "local_composition_at_spin_centroid": final_composition,
            "positive_spin_fraction_by_component": dict(sorted(component_spin.items())),
            "composition_degenerate_reference": seed["composition_degenerate_reference"],
            "inputs": {
                "neutral_input": {
                    "path": str(neutral_input.resolve()),
                    "sha256": sha256_file(neutral_input),
                },
                "neutral_output": {
                    "path": str(neutral_output.resolve()),
                    "sha256": sha256_file(neutral_output),
                },
                "anion_input": {
                    "path": str(anion_input.resolve()),
                    "sha256": sha256_file(anion_input),
                },
                "anion_output": {
                    "path": str(anion_output.resolve()),
                    "sha256": sha256_file(anion_output),
                },
                "spin_cube": {"path": str(spin_path.resolve()), "sha256": sha256_file(spin_path)},
                "electron_cube": {
                    "path": str(electron_path.resolve()),
                    "sha256": sha256_file(electron_path),
                },
                "seed_metadata": seed["metadata"],
                "spec": {"path": str(Path(args.spec).resolve()), "sha256": sha256_file(args.spec)},
            },
            **localization,
        },
    )
    return 0


def _selection_fraction(record: dict[str, Any], selection_radius: float) -> float | None:
    shells = record["local_composition_at_seed"]["shells"]
    matching = [shell for shell in shells if float(shell["radius_angstrom"]) == selection_radius]
    if len(matching) != 1:
        raise ValueError(f"selection radius {selection_radius} is absent from local composition")
    value = matching[0]["amine_mole_fraction"]
    return float(value) if value is not None else None


def run_summary(args: argparse.Namespace) -> int:
    records = [_read_json(path) for path in args.records]
    expected = {
        (system, int(replica), seed_role)
        for system in args.expected_systems
        for replica in args.expected_replicas
        for seed_role in args.expected_seed_roles
    }
    observed = {
        (
            str(record.get("system_id")),
            int(record.get("replica", -1)),
            str(record.get("seed_role")),
        )
        for record in records
    }
    problems: list[str] = []
    if observed != expected:
        problems.append(
            "observed system/replica/seed-role keys differ from expected: "
            f"{len(observed)} versus {len(expected)}"
        )
    if not records or not all(record.get("ready") for record in records):
        problems.append("one or more Stage-B2C pair analyses are not ready")
    method = _read_json(args.methods)["stage_b2c_preferential_smoke"]
    tolerance_ev = float(method["energy_tie_tolerance_ev"])
    selection_radius = float(method["selection_radius_angstrom"])
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[(str(record["system_id"]), int(record["replica"]))].append(record)
    comparisons: list[dict[str, Any]] = []
    for (system_id, replica), group in sorted(grouped.items()):
        by_role = {str(record["seed_role"]): record for record in group}
        if set(by_role) != {"eda_rich", "eda_poor"}:
            problems.append(f"{system_id} r{replica} lacks the fixed rich/poor pair")
            continue
        rich = by_role["eda_rich"]
        poor = by_role["eda_poor"]
        rich_fraction = _selection_fraction(rich, selection_radius)
        poor_fraction = _selection_fraction(poor, selection_radius)
        if (
            rich_fraction is not None
            and poor_fraction is not None
            and rich_fraction < poor_fraction
        ):
            problems.append(
                f"{system_id} r{replica} rich seed has lower EDA fraction than poor seed"
            )
        rich_energy = rich["energies"]["vertical_attachment_proxy_ev"]
        poor_energy = poor["energies"]["vertical_attachment_proxy_ev"]
        delta = (
            float(rich_energy) - float(poor_energy)
            if rich_energy is not None and poor_energy is not None
            else None
        )
        composition_degenerate = bool(
            rich.get("composition_degenerate_reference")
            and poor.get("composition_degenerate_reference")
        )
        if composition_degenerate:
            preference = "composition_degenerate_reference"
        elif delta is None:
            preference = "unavailable"
        elif abs(delta) <= tolerance_ev:
            preference = "indistinguishable_at_smoke_tolerance"
        elif delta > 0:
            preference = "eda_rich"
        else:
            preference = "eda_poor"
        comparisons.append(
            {
                "system_id": system_id,
                "amine": rich.get("amine"),
                "replica": replica,
                "ready": bool(rich.get("ready") and poor.get("ready")),
                "selection_amine_mole_fraction": {
                    "eda_rich": rich_fraction,
                    "eda_poor": poor_fraction,
                },
                "vertical_attachment_proxy_ev": {
                    "eda_rich": rich_energy,
                    "eda_poor": poor_energy,
                    "rich_minus_poor": delta,
                },
                "preference_at_smoke_tolerance": preference,
                "energy_tie_tolerance_ev": tolerance_ev,
                "composition_degenerate_reference": composition_degenerate,
                "positive_spin_fraction_by_component": {
                    "eda_rich": rich.get("positive_spin_fraction_by_component", {}),
                    "eda_poor": poor.get("positive_spin_fraction_by_component", {}),
                },
            }
        )
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "campaign": args.campaign,
            "kind": "stage_b2c_preferential_solvation_smoke",
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "problems": problems,
            "record_count": len(records),
            "expected_record_count": len(expected),
            "comparisons": comparisons,
            "records": records,
            "interpretation": method["interpretation"],
            "method_limitations": [
                "fixed nuclei test pre-existing solvent fluctuations only; "
                "no electron-induced nuclear reorganization is sampled",
                "semilocal PBE may spuriously delocalize an excess electron",
                "the charged periodic anion uses a compensating background and is not "
                "finite-size corrected",
                "ghost centers supply basis functions and initial-condition diversity; "
                "they do not restrain the electron",
                "the uncorrected vertical attachment proxy is comparable only within each "
                "matched neutral/anion geometry",
                "PFAS and Li are absent, so no Li-ionization or PFAS capture kinetics "
                "conclusion is permitted",
            ],
            "methods": {
                "path": str(Path(args.methods).resolve()),
                "sha256": sha256_file(args.methods),
            },
        },
    )
    return 0


def run_gate(args: argparse.Namespace) -> int:
    summary_path = Path(args.summary)
    summary = _read_json(summary_path)
    if not summary.get("ready"):
        raise ValueError(f"Stage-B2C summary is not ready: {summary_path}")
    digest = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(f"sha256 {digest}  {summary_path.name}\n", encoding="utf-8")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--campaign", required=True)
    prepare.add_argument("--source-campaign", required=True)
    prepare.add_argument("--candidate-manifest", required=True)
    prepare.add_argument("--candidate-summary", required=True)
    prepare.add_argument("--candidate-gate", required=True)
    prepare.add_argument("--spec", required=True)
    prepare.add_argument("--methods", required=True)
    prepare.add_argument("--systems", required=True)
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
    analyze.add_argument("--campaign", required=True)
    analyze.add_argument("--system", required=True)
    analyze.add_argument("--replica", required=True, type=int)
    analyze.add_argument("--seed-role", required=True)
    analyze.add_argument("--manifest", required=True)
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
    summary = subparsers.add_parser("summary")
    summary.add_argument("--campaign", required=True)
    summary.add_argument("--methods", required=True)
    summary.add_argument("--records", nargs="+", required=True)
    summary.add_argument("--expected-systems", nargs="+", required=True)
    summary.add_argument("--expected-replicas", nargs="+", required=True, type=int)
    summary.add_argument("--expected-seed-roles", nargs="+", required=True)
    summary.add_argument("--output", required=True)
    summary.set_defaults(func=run_summary)
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
