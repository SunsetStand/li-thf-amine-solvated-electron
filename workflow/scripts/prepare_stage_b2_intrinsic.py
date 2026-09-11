#!/usr/bin/env python3
"""Prepare, render, analyze, and gate Li-free Stage-B2 excess-electron smoke jobs."""

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
from solvelec.provenance import sha256_file
from solvelec.rendering import render_stage_b2_intrinsic_cp2k

SCIENTIFIC_STATUS = "NUMERICAL_VERTICAL_ELECTRON_ONLY_SMOKE"


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
        raise ValueError("Stage-B2 requires the accepted Stage-B candidate summary and gate")
    if not _read_json(summary).get("ready"):
        raise ValueError("accepted Stage-B candidate summary is not ready")
    digest = hashlib.sha256(summary.read_bytes()).hexdigest()
    if gate.read_text(encoding="utf-8").split() != ["sha256", digest, summary.name]:
        raise ValueError("Stage-B candidate gate does not match its summary")


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
        raise ValueError("Stage-B2 spec has no positive solvent component count")
    return counts


def _component_definitions(systems: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        "thf": systems["thf"],
        **{str(name): value for name, value in systems.get("amines", {}).items()},
    }


def _ready_seed(manifest: dict[str, Any], seed_id: str) -> dict[str, Any]:
    records = [
        record
        for record in manifest.get("seeds", [])
        if record.get("seed_id") == seed_id and record.get("ready")
    ]
    if not manifest.get("ready") or len(records) != 1:
        raise ValueError(f"Stage-B2 seed {seed_id!r} is not uniquely ready")
    return records[0]


def run_prepare(args: argparse.Namespace) -> int:
    _validate_checksum_gate(args.candidate_summary, args.candidate_gate)
    candidate_manifest_path = Path(args.candidate_manifest).resolve()
    spec_path = Path(args.spec).resolve()
    methods_path = Path(args.methods).resolve()
    candidate_manifest = _read_json(candidate_manifest_path)
    spec = _read_json(spec_path)
    methods = _read_json(methods_path)
    settings = methods["stage_b2_intrinsic_smoke"]
    if settings.get("scientific_status") != SCIENTIFIC_STATUS:
        raise ValueError("Stage-B2 intrinsic method is not explicitly smoke-only")
    if not candidate_manifest.get("ready"):
        raise ValueError("Stage-B candidate manifest is not ready")
    if str(candidate_manifest.get("system_id")) != str(spec.get("system_id")) or int(
        candidate_manifest.get("replica", -1)
    ) != int(spec.get("replica", -2)):
        raise ValueError("Stage-B candidate manifest and solvent spec differ")

    source_candidates = [
        record for record in candidate_manifest.get("candidates", []) if record.get("ready")
    ]
    if not source_candidates:
        raise ValueError("Stage-B candidate manifest contains no ready structures")
    source_structure = source_candidates[0]["structure"]
    source_xyz = Path(source_structure["xyz"]["path"])
    source_cell = Path(source_structure["cell"]["path"])
    if sha256_file(source_xyz) != source_structure["xyz"]["sha256"]:
        raise ValueError("Stage-B candidate XYZ checksum mismatch")
    if sha256_file(source_cell) != source_structure["cell"]["sha256"]:
        raise ValueError("Stage-B candidate cell checksum mismatch")
    elements, positions, _comment = read_xyz(source_xyz)
    li_indices = [index for index, element in enumerate(elements) if element.upper() == "LI"]
    ghost_indices = [index for index, element in enumerate(elements) if element.upper() == "GH"]
    if len(li_indices) != 1 or len(ghost_indices) != 1:
        raise ValueError("Stage-B2 source candidate must contain exactly one Li and one ghost")
    solvent_indices = [
        index for index, element in enumerate(elements) if element.upper() not in {"LI", "GH"}
    ]
    solvent_elements = [elements[index] for index in solvent_indices]
    solvent_positions = positions[solvent_indices]

    seed_count = int(settings["void_seed_count"])
    sites = list(candidate_manifest.get("void_search", {}).get("sites", []))
    if len(sites) < seed_count:
        raise ValueError(f"candidate manifest has {len(sites)} voids; {seed_count} are required")
    output_dir = Path(args.output_dir).resolve()
    records: list[dict[str, Any]] = []
    for seed_index, site in enumerate(sites[:seed_count], start=1):
        rank = int(site["rank"])
        seed_id = f"void_{seed_index:02d}"
        seed_dir = output_dir / seed_id
        coordinates = seed_dir / "coordinates.xyz"
        cell = seed_dir / "cell.inc"
        metadata = seed_dir / "metadata.json"
        seed_position = np.asarray(site["cartesian_angstrom"], dtype=float)
        seed_elements = [*solvent_elements, "Gh"]
        seed_positions = np.vstack([solvent_positions, seed_position])
        write_xyz(
            coordinates,
            seed_elements,
            seed_positions,
            (
                f"stage_b2_intrinsic_seed={seed_id} source={candidate_manifest['snapshot_id']} "
                "charge=-1 Li_absent PFAS_absent Gh=basis_seed_not_a_restraint"
            ),
        )
        cell.parent.mkdir(parents=True, exist_ok=True)
        cell.write_text(source_cell.read_text(encoding="utf-8"), encoding="utf-8")
        record = {
            "schema_version": 1,
            "ready": True,
            "scientific_status": SCIENTIFIC_STATUS,
            "system_id": spec["system_id"],
            "amine": spec.get("amine"),
            "replica": int(spec["replica"]),
            "snapshot_id": candidate_manifest["snapshot_id"],
            "seed_id": seed_id,
            "void_rank": rank,
            "cavity_basis_site": site,
            "electronic_state": {
                "charge": int(settings["charge"]),
                "multiplicity": int(settings["multiplicity"]),
                "li_atom_present": False,
                "pfas_present": False,
                "localization_constraint_present": False,
                "fixed_nuclei": True,
            },
            "experiment_context": settings["experiment_context"],
            "interpretation": settings["interpretation"],
            "structure": {
                "solvent_atom_count": len(solvent_elements),
                "atom_count_including_ghost": len(seed_elements),
                "ghost_basis_center_count": 1,
                "xyz": {"path": str(coordinates), "sha256": sha256_file(coordinates)},
                "cell": {"path": str(cell), "sha256": sha256_file(cell)},
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
            },
        }
        _write_json(metadata, record)
        record["metadata"] = {"path": str(metadata), "sha256": sha256_file(metadata)}
        records.append(record)

    output = {
        "schema_version": 1,
        "ready": len(records) == seed_count,
        "scientific_status": SCIENTIFIC_STATUS,
        "system_id": spec["system_id"],
        "amine": spec.get("amine"),
        "replica": int(spec["replica"]),
        "snapshot_id": candidate_manifest["snapshot_id"],
        "seed_count": len(records),
        "seed_ids": [record["seed_id"] for record in records],
        "experiment_context": settings["experiment_context"],
        "interpretation": settings["interpretation"],
        "seeds": records,
    }
    _write_json(args.output, output)
    return 0


def run_render(args: argparse.Namespace) -> int:
    manifest = _read_json(args.manifest)
    methods = _read_json(args.methods)
    seed = _ready_seed(manifest, args.seed)
    structure = seed["structure"]
    coordinates = Path(structure["xyz"]["path"])
    cell = Path(structure["cell"]["path"])
    if sha256_file(coordinates) != structure["xyz"]["sha256"]:
        raise ValueError("Stage-B2 coordinate checksum mismatch")
    if sha256_file(cell) != structure["cell"]["sha256"]:
        raise ValueError("Stage-B2 cell checksum mismatch")
    elements, _positions, _comment = read_xyz(coordinates)
    if any(element.upper() == "LI" for element in elements):
        raise ValueError("Stage-B2 intrinsic coordinates unexpectedly contain Li")
    render_stage_b2_intrinsic_cp2k(
        args.template,
        args.output,
        project=args.project,
        coordinates_path=coordinates,
        cell_path=cell,
        method=methods["stage_b2_intrinsic_smoke"],
    )
    return 0


def run_analyze(args: argparse.Namespace) -> int:
    manifest = _read_json(args.manifest)
    seed = _ready_seed(manifest, args.seed)
    spec = _read_json(args.spec)
    methods = _read_json(args.methods)
    systems = _read_json(args.systems)
    method = methods["stage_b2_intrinsic_smoke"]
    structure = seed["structure"]
    coordinates = Path(structure["xyz"]["path"])
    cell_path = Path(structure["cell"]["path"])
    for path, expected in (
        (coordinates, structure["xyz"]["sha256"]),
        (cell_path, structure["cell"]["sha256"]),
    ):
        if sha256_file(path) != expected:
            raise ValueError(f"Stage-B2 immutable structure checksum mismatch: {path}")
    if str(seed["system_id"]) != args.system or int(seed["replica"]) != int(args.replica):
        raise ValueError("Stage-B2 seed does not match requested system/replica")
    if str(spec.get("system_id")) != args.system or int(spec.get("replica", -1)) != int(
        args.replica
    ):
        raise ValueError("Stage-B2 solvent spec does not match requested system/replica")

    elements, positions, _comment = read_xyz(coordinates)
    if any(element.upper() == "LI" for element in elements):
        raise ValueError("Stage-B2 intrinsic analysis refuses a structure containing Li")
    if sum(element.upper() == "GH" for element in elements) != 1:
        raise ValueError("Stage-B2 intrinsic structure must contain exactly one ghost center")
    cell = read_cp2k_cell(cell_path)
    localization_settings = {
        **methods["stage_b_localization"],
        "cavity_probe_radius_angstrom": float(method["cavity_probe_radius_angstrom"]),
    }
    topology = infer_molecular_topology(
        elements,
        positions,
        cell,
        _component_counts_from_spec(spec),
        _component_definitions(systems),
        bond_scale=float(localization_settings["covalent_bond_scale"]),
    )
    cp2k_input_path = Path(args.cp2k_input)
    cp2k_output_path = Path(args.cp2k_output)
    spin_path = Path(args.spin_cube)
    electron_path = Path(args.electron_cube)
    input_text = cp2k_input_path.read_text(encoding="utf-8")
    upper_input = input_text.upper()
    problems = [*topology.problems]
    if not re.search(r"^\s*CHARGE\s+-1\s*$", input_text, flags=re.MULTILINE):
        problems.append("CP2K input is not an explicit charge -1 state")
    if not re.search(r"^\s*MULTIPLICITY\s+2\s*$", input_text, flags=re.MULTILINE):
        problems.append("CP2K input is not an explicit doublet")
    if "&CDFT" in upper_input or "&CONSTRAINT" in upper_input:
        problems.append("CP2K input contains a forbidden localization constraint")
    if "&KIND LI" in upper_input:
        problems.append("CP2K input contains a forbidden Li kind")
    output_text = cp2k_output_path.read_text(encoding="utf-8", errors="replace")
    result = parse_cp2k_text(output_text)
    problems.extend(result.problems)
    spin_cube = read_cube(spin_path)
    electron_cube = read_cube(electron_path)
    tolerance = float(localization_settings["cube_grid_tolerance_angstrom"])
    problems.extend(
        cube_compatibility_problems(spin_cube, electron_cube, tolerance_angstrom=tolerance)
    )
    problems.extend(
        candidate_cube_problems(spin_cube, elements, positions, cell, tolerance_angstrom=tolerance)
    )
    cavity_center = np.asarray(seed["cavity_basis_site"]["cartesian_angstrom"], dtype=float)
    localization = state_localization_record(
        spin_cube,
        elements,
        positions,
        cell,
        topology,
        cavity_center,
        localization_settings,
    )
    signed = float(localization["spin_density"]["signed_integral"])
    lower = float(localization_settings["signed_spin_integral_min"])
    upper = float(localization_settings["signed_spin_integral_max"])
    if not lower <= signed <= upper:
        problems.append(f"signed spin integral {signed:.6g} is outside [{lower:.6g}, {upper:.6g}]")

    _write_json(
        args.output,
        {
            "schema_version": 1,
            "campaign": args.campaign,
            "kind": "stage_b2_intrinsic_state",
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "system_id": args.system,
            "amine": spec.get("amine"),
            "replica": int(args.replica),
            "seed_id": args.seed,
            "void_rank": int(seed["void_rank"]),
            "energy_hartree": result.energy_hartree,
            "normal_termination": result.normal_termination,
            "scf_converged": result.converged,
            "problems": problems,
            "checks": {
                "li_absent": True,
                "pfas_absent": True,
                "charge_minus_one_doublet": not any(
                    "explicit charge" in problem or "explicit doublet" in problem
                    for problem in problems
                ),
                "constraint_release": not any("constraint" in problem for problem in problems),
                "cube_grid_and_atoms": not any(problem.startswith("cube") for problem in problems),
                "molecular_topology": not topology.problems,
                "signed_spin_integral": not any(
                    problem.startswith("signed spin integral") for problem in problems
                ),
            },
            "experiment_context": method["experiment_context"],
            "interpretation": method["interpretation"],
            "inputs": {
                "cp2k_input": {
                    "path": str(cp2k_input_path.resolve()),
                    "sha256": sha256_file(cp2k_input_path),
                },
                "cp2k_output": {
                    "path": str(cp2k_output_path.resolve()),
                    "sha256": sha256_file(cp2k_output_path),
                },
                "spin_cube": {"path": str(spin_path.resolve()), "sha256": sha256_file(spin_path)},
                "electron_cube": {
                    "path": str(electron_path.resolve()),
                    "sha256": sha256_file(electron_path),
                },
                "seed_metadata": seed["metadata"],
                "coordinates": {
                    "path": str(coordinates.resolve()),
                    "sha256": sha256_file(coordinates),
                },
                "cell": {"path": str(cell_path.resolve()), "sha256": sha256_file(cell_path)},
                "spec": {"path": str(Path(args.spec).resolve()), "sha256": sha256_file(args.spec)},
                "methods": {
                    "path": str(Path(args.methods).resolve()),
                    "sha256": sha256_file(args.methods),
                },
            },
            **localization,
        },
    )
    return 0


def _proxy_signature(record: dict[str, Any]) -> str:
    flags = record.get("localization_proxy_flags", {})
    active = sorted(name for name, value in flags.items() if value and name != "li_centered")
    return "+".join(active) if active else "none"


def run_summary(args: argparse.Namespace) -> int:
    records = [_read_json(path) for path in args.records]
    expected = {
        (system, int(replica), seed)
        for system in args.expected_systems
        for replica in args.expected_replicas
        for seed in args.expected_seeds
    }
    observed = {
        (str(record.get("system_id")), int(record.get("replica", -1)), str(record.get("seed_id")))
        for record in records
    }
    problems: list[str] = []
    if observed != expected:
        problems.append(
            "observed system/replica/seed keys differ from expected: "
            f"{len(observed)} versus {len(expected)}"
        )
    if not records or not all(record.get("ready") for record in records):
        problems.append("one or more Stage-B2 intrinsic analyses are not ready")
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[(str(record["system_id"]), int(record["replica"]))].append(record)
    systems: list[dict[str, Any]] = []
    for (system_id, replica), group in sorted(grouped.items()):
        energies = [
            float(record["energy_hartree"])
            for record in group
            if record["energy_hartree"] is not None
        ]
        signatures = sorted({_proxy_signature(record) for record in group})
        systems.append(
            {
                "system_id": system_id,
                "amine": group[0].get("amine"),
                "replica": replica,
                "ready": all(record.get("ready") for record in group),
                "seed_ids": sorted(str(record["seed_id"]) for record in group),
                "energy_spread_hartree": max(energies) - min(energies) if energies else None,
                "energy_spread_ev": (
                    (max(energies) - min(energies)) * HARTREE_TO_EV if energies else None
                ),
                "proxy_signatures": signatures,
                "seed_sensitive_proxy": len(signatures) > 1,
            }
        )
    methods = _read_json(args.methods)
    method = methods["stage_b2_intrinsic_smoke"]
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "campaign": args.campaign,
            "kind": "stage_b2_intrinsic_smoke",
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "problems": problems,
            "record_count": len(records),
            "expected_record_count": len(expected),
            "experiment_context": method["experiment_context"],
            "interpretation": method["interpretation"],
            "method_limitations": [
                "fixed nuclei: this is a vertical injection state, not relaxed electron solvation",
                "semilocal PBE may spuriously delocalize an excess electron",
                "one charged periodic cell size is not a finite-size convergence test",
                "ghost centers are basis seeds; their probe spheres are not "
                "observed physical cavities",
                "PFAS is absent, so no capture or reaction-rate conclusion is permitted",
            ],
            "systems": systems,
            "records": records,
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
        raise ValueError(f"Stage-B2 intrinsic summary is not ready: {summary_path}")
    digest = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(f"sha256 {digest}  {summary_path.name}\n", encoding="utf-8")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--candidate-manifest", required=True)
    prepare.add_argument("--candidate-summary", required=True)
    prepare.add_argument("--candidate-gate", required=True)
    prepare.add_argument("--spec", required=True)
    prepare.add_argument("--methods", required=True)
    prepare.add_argument("--output-dir", required=True)
    prepare.add_argument("--output", required=True)
    prepare.set_defaults(func=run_prepare)
    render = subparsers.add_parser("render")
    render.add_argument("--manifest", required=True)
    render.add_argument("--methods", required=True)
    render.add_argument("--template", required=True)
    render.add_argument("--seed", required=True)
    render.add_argument("--project", required=True)
    render.add_argument("--output", required=True)
    render.set_defaults(func=run_render)
    analyze = subparsers.add_parser("analyze")
    analyze.add_argument("--campaign", required=True)
    analyze.add_argument("--system", required=True)
    analyze.add_argument("--replica", required=True, type=int)
    analyze.add_argument("--seed", required=True)
    analyze.add_argument("--manifest", required=True)
    analyze.add_argument("--spec", required=True)
    analyze.add_argument("--methods", required=True)
    analyze.add_argument("--systems", required=True)
    analyze.add_argument("--cp2k-input", required=True)
    analyze.add_argument("--cp2k-output", required=True)
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
    summary.add_argument("--expected-seeds", nargs="+", required=True)
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
