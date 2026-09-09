#!/usr/bin/env python3
"""Prepare and gate generic Li/cavity candidates from immutable Stage-A snapshots."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from solvelec.candidates import (
    heavy_atom_geometry,
    ranked_void_sites,
    read_xyz,
    select_candidate_pairs,
    write_xyz,
)
from solvelec.parsers import HARTREE_TO_EV, evaluate_cp2k_cdft_constraint, parse_cp2k_text
from solvelec.provenance import sha256_file
from solvelec.rendering import render_stage_b_cp2k


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


def _validate_handoff(
    spec: dict[str, Any], metadata: dict[str, Any], xyz_path: Path, cell_path: Path
) -> np.ndarray:
    if not metadata.get("ready") or metadata.get("li_atom_present") is not False:
        raise ValueError("Stage-A snapshot must be ready and explicitly Li-free")
    if str(metadata.get("system_id")) != str(spec.get("system_id")):
        raise ValueError("spec and snapshot system ids differ")
    if int(metadata.get("replica", -1)) != int(spec.get("replica", -2)):
        raise ValueError("spec and snapshot replicas differ")
    structure = metadata.get("structure", {})
    expected_xyz = structure.get("xyz", {}).get("sha256")
    expected_cell = structure.get("cell", {}).get("sha256")
    if expected_xyz != sha256_file(xyz_path) or expected_cell != sha256_file(cell_path):
        raise ValueError("Stage-A snapshot checksum mismatch")
    matrix = np.asarray(structure.get("cell_vectors_angstrom"), dtype=float)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError("snapshot metadata does not contain a finite 3x3 cell")
    if float(np.linalg.det(matrix)) <= 0:
        raise ValueError("snapshot cell must have positive volume")
    return matrix


def run_prepare(args: argparse.Namespace) -> int:
    spec_path = Path(args.spec).resolve()
    snapshot_metadata_path = Path(args.snapshot_metadata).resolve()
    xyz_path = Path(args.xyz).resolve()
    cell_path = Path(args.cell).resolve()
    methods_path = Path(args.methods).resolve()
    output_dir = Path(args.output_dir).resolve()
    spec = _read_json(spec_path)
    metadata = _read_json(snapshot_metadata_path)
    methods = _read_json(methods_path)
    settings = methods["stage_b"]
    matrix = _validate_handoff(spec, metadata, xyz_path, cell_path)
    elements, positions, _ = read_xyz(xyz_path)
    if any(element.upper() == "LI" for element in elements):
        raise ValueError("Stage-A source XYZ unexpectedly contains Li")
    if len(elements) != int(metadata["structure"]["atom_count"]):
        raise ValueError("Stage-A XYZ atom count disagrees with metadata")
    heavy_positions, radii = heavy_atom_geometry(elements, positions)
    sites = ranked_void_sites(
        heavy_positions,
        radii,
        matrix,
        count=int(settings["candidate_site_count"]),
        points_per_axis=int(settings["void_grid_points_per_axis"]),
        refinement_levels=int(settings["void_refinement_levels"]),
        minimum_separation_angstrom=float(settings["minimum_site_separation_angstrom"]),
        minimum_clearance_angstrom=float(settings["minimum_surface_clearance_angstrom"]),
    )
    candidates = select_candidate_pairs(sites, matrix, settings["candidate_pairs"])
    records: list[dict[str, Any]] = []
    for candidate in candidates:
        candidate_dir = output_dir / str(candidate["id"])
        candidate_xyz = candidate_dir / "coordinates.xyz"
        candidate_cell = candidate_dir / "cell.inc"
        candidate_metadata = candidate_dir / "metadata.json"
        li_position = np.asarray(candidate["li_site"]["cartesian_angstrom"], dtype=float)
        ghost_position = np.asarray(
            candidate["cavity_basis_site"]["cartesian_angstrom"], dtype=float
        )
        candidate_elements = ["LI", *elements, "Gh"]
        candidate_positions = np.vstack([li_position, positions, ghost_position])
        write_xyz(
            candidate_xyz,
            candidate_elements,
            candidate_positions,
            (
                f"stage_b_candidate={candidate['id']} source={metadata['snapshot_id']} "
                "Gh=cavity_basis_seed_not_an_atom"
            ),
        )
        candidate_cell.parent.mkdir(parents=True, exist_ok=True)
        candidate_cell.write_text(cell_path.read_text(encoding="utf-8"), encoding="utf-8")
        record = {
            "schema_version": 1,
            "ready": True,
            "scientific_status": "GEOMETRIC_SEED_NOT_LOCALIZATION_EVIDENCE",
            "system_id": spec["system_id"],
            "amine": spec.get("amine"),
            "replica": int(spec["replica"]),
            "snapshot_id": metadata["snapshot_id"],
            "candidate_id": candidate["id"],
            "target_li_cavity_distance_angstrom": candidate["target_li_cavity_distance_angstrom"],
            "tolerance_angstrom": candidate["tolerance_angstrom"],
            "achieved_li_cavity_distance_angstrom": candidate[
                "achieved_li_cavity_distance_angstrom"
            ],
            "li_site": candidate["li_site"],
            "cavity_basis_site": candidate["cavity_basis_site"],
            "structure": {
                "atom_count_including_li_and_ghost": len(candidate_elements),
                "solvent_atom_count": len(elements),
                "li_atom_index_cp2k": 1,
                "ghost_basis_center_count": 1,
                "xyz": {
                    "path": str(candidate_xyz),
                    "sha256": sha256_file(candidate_xyz),
                },
                "cell": {
                    "path": str(candidate_cell),
                    "sha256": sha256_file(candidate_cell),
                },
            },
            "source": {
                "spec": {"path": str(spec_path), "sha256": sha256_file(spec_path)},
                "methods": {
                    "path": str(methods_path),
                    "sha256": sha256_file(methods_path),
                },
                "snapshot_metadata": {
                    "path": str(snapshot_metadata_path),
                    "sha256": sha256_file(snapshot_metadata_path),
                },
                "snapshot_xyz_sha256": sha256_file(xyz_path),
                "snapshot_cell_sha256": sha256_file(cell_path),
            },
        }
        _write_json(candidate_metadata, record)
        record["metadata"] = {
            "path": str(candidate_metadata),
            "sha256": sha256_file(candidate_metadata),
        }
        records.append(record)
    manifest = {
        "schema_version": 1,
        "ready": len(records) == len(settings["candidate_pairs"]),
        "scientific_status": "CANDIDATE_BANK_ONLY",
        "system_id": spec["system_id"],
        "amine": spec.get("amine"),
        "replica": int(spec["replica"]),
        "snapshot_id": metadata["snapshot_id"],
        "candidate_count": len(records),
        "candidate_ids": [record["candidate_id"] for record in records],
        "void_search": {
            "site_count": len(sites),
            "sites": sites,
            "heavy_atom_count": len(heavy_positions),
        },
        "configuration": settings,
        "methods": {"path": str(methods_path), "sha256": sha256_file(methods_path)},
        "candidates": records,
    }
    _write_json(args.output, manifest)
    return 0


def run_render(args: argparse.Namespace) -> int:
    manifest = _read_json(args.manifest)
    methods = _read_json(args.methods)
    matches = [
        record
        for record in manifest.get("candidates", [])
        if record.get("candidate_id") == args.candidate
    ]
    if not manifest.get("ready") or len(matches) != 1:
        raise ValueError(f"candidate {args.candidate!r} is not uniquely ready")
    candidate = matches[0]
    structure = candidate["structure"]
    xyz_path = Path(structure["xyz"]["path"])
    cell_path = Path(structure["cell"]["path"])
    if sha256_file(xyz_path) != structure["xyz"]["sha256"]:
        raise ValueError("candidate XYZ checksum mismatch")
    if sha256_file(cell_path) != structure["cell"]["sha256"]:
        raise ValueError("candidate cell checksum mismatch")
    render_stage_b_cp2k(
        args.template,
        args.output,
        project=args.project,
        coordinates_path=xyz_path,
        cell_path=cell_path,
        method=methods["stage_b_smoke"],
        li_atom_index=int(structure["li_atom_index_cp2k"]),
        target_electrons=args.target_electrons,
    )
    return 0


def summarize_manifests(records: list[dict[str, Any]], campaign: str) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[str(record.get("system_id"))].append(record)
    systems: dict[str, Any] = {}
    for system_id, group in sorted(grouped.items()):
        systems[system_id] = {
            "amine": group[0].get("amine"),
            "replicas": sorted(int(record["replica"]) for record in group),
            "ready": all(bool(record.get("ready")) for record in group),
            "candidate_ids": sorted(
                {candidate for record in group for candidate in record.get("candidate_ids", [])}
            ),
        }
    keys = [(str(record.get("system_id")), int(record.get("replica", -1))) for record in records]
    candidate_sets = [tuple(record.get("candidate_ids", [])) for record in records]
    unique_records = len(keys) == len(set(keys))
    consistent_candidates = bool(candidate_sets) and len(set(candidate_sets)) == 1
    return {
        "schema_version": 1,
        "campaign": campaign,
        "kind": "stage_b_candidates",
        "ready": (
            bool(records)
            and unique_records
            and consistent_candidates
            and all(bool(record.get("ready")) for record in records)
        ),
        "record_count": len(records),
        "unique_system_replica_records": unique_records,
        "consistent_candidate_ids": consistent_candidates,
        "systems": systems,
        "scientific_status": "CANDIDATE_BANK_ONLY",
    }


def run_summary(args: argparse.Namespace) -> int:
    records = [_read_json(path) for path in args.inputs]
    _write_json(args.output, summarize_manifests(records, args.campaign))
    return 0


def run_smoke_summary(args: argparse.Namespace) -> int:
    if not (len(args.outputs) == len(args.manifests) == len(args.cp2k_inputs)):
        raise ValueError("smoke outputs, inputs, and manifests must have the same cardinality")
    methods = _read_json(args.methods)
    smoke_method = methods["stage_b_smoke"]
    cdft_tolerance = float(smoke_method["cdft_eps_scf"])
    cdft_target = float(smoke_method["li_target_valence_electrons"])
    records: list[dict[str, Any]] = []
    for output_path, input_path, manifest_path in zip(
        args.outputs, args.cp2k_inputs, args.manifests, strict=True
    ):
        output = Path(output_path)
        cp2k_input = Path(input_path)
        output_text = output.read_text(encoding="utf-8", errors="replace")
        result = parse_cp2k_text(output_text)
        cdft = evaluate_cp2k_cdft_constraint(
            output_text,
            expected_target_electrons=cdft_target,
            tolerance_electrons=cdft_tolerance,
        )
        manifest = _read_json(manifest_path)
        matches = [
            record
            for record in manifest.get("candidates", [])
            if record.get("candidate_id") == args.candidate and record.get("ready")
        ]
        if len(matches) != 1:
            raise ValueError(
                f"manifest {manifest_path} lacks one ready {args.candidate!r} candidate"
            )
        records.append(
            {
                "system_id": manifest["system_id"],
                "amine": manifest.get("amine"),
                "replica": int(manifest["replica"]),
                "candidate_id": args.candidate,
                "converged": result.converged and cdft.converged,
                "normal_termination": result.normal_termination,
                "energy_hartree": result.energy_hartree,
                "problems": [*result.problems, *cdft.problems],
                "cdft_constraint_gate": cdft.as_dict(),
                "input": {
                    "path": str(cp2k_input.resolve()),
                    "sha256": sha256_file(cp2k_input),
                },
                "output": {"path": str(output.resolve()), "sha256": sha256_file(output)},
            }
        )
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "campaign": args.campaign,
            "kind": "stage_b_cp2k_smoke",
            "ready": bool(records) and all(record["converged"] for record in records),
            "scientific_status": "NUMERICAL_SMOKE_ONLY_NOT_A_LOCALIZATION_RESULT",
            "records": records,
        },
    )
    return 0


def run_mechanism_summary(args: argparse.Namespace) -> int:
    cardinalities = {
        len(args.outputs),
        len(args.cp2k_inputs),
        len(args.manifests),
        len(args.states),
        len(args.targets),
    }
    if len(cardinalities) != 1:
        raise ValueError(
            "mechanism outputs, inputs, manifests, states, and targets must have "
            "the same cardinality"
        )
    methods = _read_json(args.methods)
    smoke_method = methods["stage_b_smoke"]
    tolerance = float(smoke_method["cdft_eps_scf"])
    state_configuration = {str(record["id"]): record for record in smoke_method["mechanism_states"]}
    configured_states = {
        state_id: float(record["li_target_valence_electrons"])
        for state_id, record in state_configuration.items()
    }
    expected_state_ids = {"li0_diabatic", "li_plus_e_diabatic"}
    if set(configured_states) != expected_state_ids:
        raise ValueError("Stage-B mechanism smoke requires the configured Li0 and Li+e states")

    states: list[dict[str, Any]] = []
    for output_path, input_path, manifest_path, state_id, target_value in zip(
        args.outputs,
        args.cp2k_inputs,
        args.manifests,
        args.states,
        args.targets,
        strict=True,
    ):
        target = float(target_value)
        if state_id not in configured_states or target != configured_states[state_id]:
            raise ValueError(f"state {state_id!r} does not match the configured cDFT target")
        output = Path(output_path)
        cp2k_input = Path(input_path)
        output_text = output.read_text(encoding="utf-8", errors="replace")
        result = parse_cp2k_text(output_text)
        cdft = evaluate_cp2k_cdft_constraint(
            output_text,
            expected_target_electrons=target,
            tolerance_electrons=tolerance,
        )
        manifest = _read_json(manifest_path)
        matches = [
            record
            for record in manifest.get("candidates", [])
            if record.get("candidate_id") == args.candidate and record.get("ready")
        ]
        if len(matches) != 1:
            raise ValueError(
                f"manifest {manifest_path} lacks one ready {args.candidate!r} candidate"
            )
        states.append(
            {
                "system_id": manifest["system_id"],
                "amine": manifest.get("amine"),
                "replica": int(manifest["replica"]),
                "candidate_id": args.candidate,
                "state_id": state_id,
                "li_target_valence_electrons": target,
                "interpretation": str(state_configuration[state_id]["interpretation"]),
                "converged": result.converged and cdft.converged,
                "normal_termination": result.normal_termination,
                "energy_hartree": result.energy_hartree,
                "problems": [*result.problems, *cdft.problems],
                "cdft_constraint_gate": cdft.as_dict(),
                "input": {
                    "path": str(cp2k_input.resolve()),
                    "sha256": sha256_file(cp2k_input),
                },
                "output": {"path": str(output.resolve()), "sha256": sha256_file(output)},
            }
        )

    grouped: dict[tuple[str, int, str], list[dict[str, Any]]] = defaultdict(list)
    for state in states:
        grouped[
            (str(state["system_id"]), int(state["replica"]), str(state["candidate_id"]))
        ].append(state)
    records: list[dict[str, Any]] = []
    for (system_id, replica, candidate_id), group in sorted(grouped.items()):
        by_state = {str(state["state_id"]): state for state in group}
        complete_pair = len(group) == 2 and set(by_state) == expected_state_ids
        li0 = by_state.get("li0_diabatic")
        separated = by_state.get("li_plus_e_diabatic")
        energies_present = (
            li0 is not None
            and separated is not None
            and li0["energy_hartree"] is not None
            and separated["energy_hartree"] is not None
        )
        delta_hartree = None
        if energies_present:
            delta_hartree = float(separated["energy_hartree"]) - float(li0["energy_hartree"])
        pair_ready = complete_pair and all(bool(state["converged"]) for state in group)
        records.append(
            {
                "system_id": system_id,
                "amine": group[0].get("amine"),
                "replica": replica,
                "candidate_id": candidate_id,
                "ready": pair_ready,
                "complete_pair": complete_pair,
                "states": sorted(group, key=lambda record: str(record["state_id"])),
                "delta_energy_hartree_li_plus_e_minus_li0": delta_hartree,
                "fixed_geometry_diabatic_gap_ev": (
                    delta_hartree * HARTREE_TO_EV if delta_hartree is not None else None
                ),
            }
        )
    expected_record_count = len(states) // 2
    ready = (
        bool(records)
        and len(states) % 2 == 0
        and len(records) == expected_record_count
        and all(record["ready"] for record in records)
    )
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "campaign": args.campaign,
            "kind": "stage_b_mechanism_smoke",
            "ready": ready,
            "scientific_status": (
                "NUMERICAL_MECHANISM_SMOKE_ONLY_NOT_A_STABILITY_OR_LOCALIZATION_RESULT"
            ),
            "energy_definition": "E(li_plus_e_diabatic) - E(li0_diabatic)",
            "methods": {
                "path": str(Path(args.methods).resolve()),
                "sha256": sha256_file(args.methods),
            },
            "cdft_tolerance_electrons": tolerance,
            "interpretation": (
                "Fixed-geometry, same-electron-count paired cDFT execution check. "
                "Neither the sign nor magnitude is a solvated-electron stability result."
            ),
            "records": records,
        },
    )
    return 0


def run_gate(args: argparse.Namespace) -> int:
    summary_path = Path(args.summary)
    summary = _read_json(summary_path)
    if not summary.get("ready"):
        raise ValueError(f"Stage-B summary is not ready: {summary_path}")
    digest = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(f"sha256 {digest}  {summary_path.name}\n", encoding="utf-8")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--spec", required=True)
    prepare.add_argument("--snapshot-metadata", required=True)
    prepare.add_argument("--xyz", required=True)
    prepare.add_argument("--cell", required=True)
    prepare.add_argument("--methods", required=True)
    prepare.add_argument("--output-dir", required=True)
    prepare.add_argument("--output", required=True)
    prepare.set_defaults(func=run_prepare)
    render = subparsers.add_parser("render")
    render.add_argument("--manifest", required=True)
    render.add_argument("--methods", required=True)
    render.add_argument("--template", required=True)
    render.add_argument("--candidate", required=True)
    render.add_argument("--project", required=True)
    render.add_argument("--target-electrons", type=float)
    render.add_argument("--output", required=True)
    render.set_defaults(func=run_render)
    summary = subparsers.add_parser("summary")
    summary.add_argument("--campaign", required=True)
    summary.add_argument("--output", required=True)
    summary.add_argument("inputs", nargs="+")
    summary.set_defaults(func=run_summary)
    smoke = subparsers.add_parser("smoke-summary")
    smoke.add_argument("--campaign", required=True)
    smoke.add_argument("--candidate", required=True)
    smoke.add_argument("--methods", required=True)
    smoke.add_argument("--output", required=True)
    smoke.add_argument("--outputs", nargs="+", required=True)
    smoke.add_argument("--cp2k-inputs", nargs="+", required=True)
    smoke.add_argument("--manifests", nargs="+", required=True)
    smoke.set_defaults(func=run_smoke_summary)
    mechanism = subparsers.add_parser("mechanism-summary")
    mechanism.add_argument("--campaign", required=True)
    mechanism.add_argument("--candidate", required=True)
    mechanism.add_argument("--methods", required=True)
    mechanism.add_argument("--output", required=True)
    mechanism.add_argument("--outputs", nargs="+", required=True)
    mechanism.add_argument("--cp2k-inputs", nargs="+", required=True)
    mechanism.add_argument("--manifests", nargs="+", required=True)
    mechanism.add_argument("--states", nargs="+", required=True)
    mechanism.add_argument("--targets", nargs="+", required=True, type=float)
    mechanism.set_defaults(func=run_mechanism_summary)
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
