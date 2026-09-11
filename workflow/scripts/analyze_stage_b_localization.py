#!/usr/bin/env python3
"""Analyze Stage-B cube files without rerunning electronic-structure jobs."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from solvelec.candidates import read_xyz
from solvelec.cube import CubeData, read_cube, write_cube
from solvelec.localization_analysis import (
    candidate_cube_problems,
    cube_compatibility_problems,
    density_difference_record,
    infer_molecular_topology,
    read_cp2k_cell,
    state_localization_record,
)
from solvelec.provenance import sha256_file


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


def _write_csv(path: str | Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)


def _component_definitions(systems: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        "thf": systems["thf"],
        **{str(name): value for name, value in systems.get("amines", {}).items()},
    }


def _component_counts_from_spec(spec: dict[str, Any]) -> dict[str, int]:
    """Read current specs and reconstruct immutable pre-schema Stage-A specs."""
    configured = spec.get("component_counts")
    if configured is not None:
        if not isinstance(configured, dict) or not configured:
            raise ValueError("component_counts must be a non-empty mapping")
        counts = {str(name): int(count) for name, count in configured.items()}
    else:
        # Pilot Stage A predates component_counts. Preserve that accepted
        # handoff and recover the equivalent counts from its original fields.
        counts = {}
        thf_count = int(spec.get("thf_count", 0))
        if thf_count > 0:
            counts["thf"] = thf_count
        amine = spec.get("amine")
        amine_count = int(spec.get("amine_count_initial", 0))
        if amine:
            if amine_count <= 0:
                raise ValueError("legacy mixed-solvent spec has no positive amine count")
            counts[str(amine)] = amine_count
        elif amine_count != 0:
            raise ValueError("legacy spec has an amine count but no amine identifier")
    if not counts or any(count <= 0 for count in counts.values()):
        raise ValueError("spec must contain at least one positive solvent component count")
    return counts


def _validated_inputs(args: argparse.Namespace) -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    list[str],
    np.ndarray,
    np.ndarray,
    Any,
]:
    metadata = _read_json(args.candidate_metadata)
    spec = _read_json(args.spec)
    methods = _read_json(args.methods)
    systems = _read_json(args.systems)
    if not metadata.get("ready"):
        raise ValueError("Stage-B candidate metadata is not ready")
    for key, expected in (
        ("system_id", args.system),
        ("replica", int(args.replica)),
        ("candidate_id", args.candidate),
    ):
        actual = metadata.get(key)
        if key == "replica":
            actual = int(actual)
        if actual != expected:
            raise ValueError(f"candidate metadata {key}={actual!r}, expected {expected!r}")
    if str(spec.get("system_id")) != args.system or int(spec.get("replica", -1)) != int(
        args.replica
    ):
        raise ValueError("Stage-B localization spec does not match requested system/replica")
    coordinates = Path(args.coordinates)
    cell_path = Path(args.cell)
    structure = metadata.get("structure", {})
    if structure.get("xyz", {}).get("sha256") != sha256_file(coordinates):
        raise ValueError("candidate XYZ checksum differs from its immutable metadata")
    if structure.get("cell", {}).get("sha256") != sha256_file(cell_path):
        raise ValueError("candidate cell checksum differs from its immutable metadata")
    elements, positions, _ = read_xyz(coordinates)
    cell = read_cp2k_cell(cell_path)
    settings = methods["stage_b_localization"]
    topology = infer_molecular_topology(
        elements,
        positions,
        cell,
        _component_counts_from_spec(spec),
        _component_definitions(systems),
        bond_scale=float(settings["covalent_bond_scale"]),
    )
    return metadata, spec, methods, elements, positions, cell, topology


def run_state(args: argparse.Namespace) -> int:
    metadata, spec, methods, elements, positions, cell, topology = _validated_inputs(args)
    settings = methods["stage_b_localization"]
    spin_path = Path(args.spin_cube)
    electron_path = Path(args.electron_cube)
    spin_cube = read_cube(spin_path)
    electron_cube = read_cube(electron_path)
    tolerance = float(settings["cube_grid_tolerance_angstrom"])
    problems = [*topology.problems]
    problems.extend(
        cube_compatibility_problems(
            spin_cube, electron_cube, tolerance_angstrom=tolerance
        )
    )
    problems.extend(
        candidate_cube_problems(
            spin_cube,
            elements,
            positions,
            cell,
            tolerance_angstrom=tolerance,
        )
    )
    cavity_center = np.asarray(metadata["cavity_basis_site"]["cartesian_angstrom"], dtype=float)
    analysis = state_localization_record(
        spin_cube,
        elements,
        positions,
        cell,
        topology,
        cavity_center,
        settings,
    )
    signed_integral = float(analysis["spin_density"]["signed_integral"])
    lower = float(settings["signed_spin_integral_min"])
    upper = float(settings["signed_spin_integral_max"])
    if not lower <= signed_integral <= upper:
        problems.append(
            f"signed spin integral {signed_integral:.6g} is outside [{lower:.6g}, {upper:.6g}]"
        )
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "campaign": args.campaign,
            "kind": "stage_b_localization_state",
            "ready": not problems,
            "scientific_status": "NUMERICAL_LOCALIZATION_DIAGNOSTIC_ONLY",
            "system_id": args.system,
            "amine": spec.get("amine"),
            "replica": int(args.replica),
            "candidate_id": args.candidate,
            "state_id": args.state,
            "problems": problems,
            "checks": {
                "signed_spin_integral": not any(
                    problem.startswith("signed spin integral") for problem in problems
                ),
                "cube_grid_and_atoms": not any(
                    problem.startswith("cube") for problem in problems
                ),
                "molecular_topology": not topology.problems,
            },
            "interpretation": (
                "Positive spin-density fractions in geometric regions are screening proxies. "
                "They are not Bader/Hirshfeld populations or proof of a stable solvated electron."
            ),
            "inputs": {
                "spin_cube": {"path": str(spin_path.resolve()), "sha256": sha256_file(spin_path)},
                "electron_cube": {
                    "path": str(electron_path.resolve()),
                    "sha256": sha256_file(electron_path),
                },
                "candidate_metadata": {
                    "path": str(Path(args.candidate_metadata).resolve()),
                    "sha256": sha256_file(args.candidate_metadata),
                },
                "coordinates": {
                    "path": str(Path(args.coordinates).resolve()),
                    "sha256": sha256_file(args.coordinates),
                },
                "cell": {"path": str(Path(args.cell).resolve()), "sha256": sha256_file(args.cell)},
                "spec": {"path": str(Path(args.spec).resolve()), "sha256": sha256_file(args.spec)},
                "methods": {
                    "path": str(Path(args.methods).resolve()),
                    "sha256": sha256_file(args.methods),
                },
            },
            **analysis,
        },
    )
    return 0


def run_pair(args: argparse.Namespace) -> int:
    metadata, spec, methods, elements, positions, _cell, topology = _validated_inputs(args)
    settings = methods["stage_b_localization"]
    state_records = [_read_json(path) for path in args.state_records]
    by_state = {str(record.get("state_id")): record for record in state_records}
    required = {"li0_diabatic", "li_plus_e_diabatic"}
    if set(by_state) != required:
        raise ValueError(f"localization pair requires states {sorted(required)}")
    li0_path = Path(args.li0_electron_cube)
    separated_path = Path(args.li_plus_e_electron_cube)
    li0_cube = read_cube(li0_path)
    separated_cube = read_cube(separated_path)
    cavity_center = np.asarray(metadata["cavity_basis_site"]["cartesian_angstrom"], dtype=float)
    difference = density_difference_record(
        li0_cube,
        separated_cube,
        elements,
        positions,
        topology,
        cavity_center,
        settings,
    )
    difference_cube_path = Path(args.difference_cube)
    if difference.get("ready"):
        write_cube(
            difference_cube_path,
            CubeData(
                comments=(
                    "Stage-B1 density redistribution: li_plus_e_diabatic - li0_diabatic",
                    "Same-electron-count diagnostic; positive means density accumulation",
                ),
                origin=li0_cube.origin,
                axes=li0_cube.axes,
                atoms=li0_cube.atoms,
                values=separated_cube.values - li0_cube.values,
                coordinate_unit=li0_cube.coordinate_unit,
            ),
        )
    problems = [*topology.problems, *difference.get("problems", [])]
    if not all(record.get("ready") for record in state_records):
        problems.append("one or more state-localization records are not ready")
    difference["problems"] = problems
    difference["ready"] = not problems
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "campaign": args.campaign,
            "kind": "stage_b_localization_pair",
            "ready": not problems,
            "scientific_status": "NUMERICAL_DENSITY_REARRANGEMENT_DIAGNOSTIC_ONLY",
            "system_id": args.system,
            "amine": spec.get("amine"),
            "replica": int(args.replica),
            "candidate_id": args.candidate,
            "interpretation": (
                "The two neutral-doublet states have the same total electron count. "
                "Their density difference describes redistribution, not addition of an electron."
            ),
            "inputs": {
                "li0_electron_cube": {
                    "path": str(li0_path.resolve()),
                    "sha256": sha256_file(li0_path),
                },
                "li_plus_e_electron_cube": {
                    "path": str(separated_path.resolve()),
                    "sha256": sha256_file(separated_path),
                },
                "state_records": [
                    {"path": str(Path(path).resolve()), "sha256": sha256_file(path)}
                    for path in args.state_records
                ],
                "density_difference_cube": (
                    {
                        "path": str(difference_cube_path.resolve()),
                        "sha256": sha256_file(difference_cube_path),
                    }
                    if difference_cube_path.is_file()
                    else None
                ),
            },
            "density_difference": difference,
        },
    )
    return 0


def _state_csv_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        spin = record["spin_density"]
        region = record["geometric_partition"]
        rows.append(
            {
                "system_id": record["system_id"],
                "replica": record["replica"],
                "candidate_id": record["candidate_id"],
                "state_id": record["state_id"],
                "ready": record["ready"],
                "signed_spin_electrons": spin["signed_integral"],
                "positive_spin_electrons": spin["electron_count"],
                "negative_spin_magnitude_electrons": spin["negative_magnitude_integral"],
                "radius_angstrom": spin["radius_angstrom"],
                "inverse_participation_ratio": spin["inverse_participation_ratio"],
                "centroid_to_li_angstrom": spin["centroid_to_li_angstrom"],
                "centroid_to_ghost_cavity_angstrom": spin[
                    "centroid_to_ghost_cavity_angstrom"
                ],
                "li_positive_spin_fraction": region["li_positive_spin_fraction"],
                "maximum_solvent_molecule_positive_spin_fraction": region[
                    "maximum_solvent_molecule_positive_spin_fraction"
                ],
                "interstitial_positive_spin_fraction": region[
                    "interstitial_positive_spin_fraction"
                ],
                "ghost_cavity_probe_positive_spin_fraction": region[
                    "ghost_cavity_probe_positive_spin_fraction"
                ],
            }
        )
    return rows


def _molecule_csv_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        for molecule in record["molecules"]:
            rows.append(
                {
                    "system_id": record["system_id"],
                    "replica": record["replica"],
                    "candidate_id": record["candidate_id"],
                    "state_id": record["state_id"],
                    "molecule_id": molecule["molecule_id"],
                    "component": molecule["component"],
                    "signed_spin_electrons": molecule["signed_spin_electrons"],
                    "positive_spin_electrons": molecule["positive_spin_electrons"],
                    "negative_spin_magnitude_electrons": molecule[
                        "negative_spin_magnitude_electrons"
                    ],
                    "positive_spin_fraction": molecule["positive_spin_fraction"],
                }
            )
    return rows


def _pair_csv_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        density = record["density_difference"]
        rows.append(
            {
                "system_id": record["system_id"],
                "replica": record["replica"],
                "candidate_id": record["candidate_id"],
                "ready": record["ready"],
                "signed_integral_electrons": density.get("signed_integral_electrons"),
                "accumulation_integral_electrons": density.get(
                    "accumulation_integral_electrons"
                ),
                "depletion_magnitude_integral_electrons": density.get(
                    "depletion_magnitude_integral_electrons"
                ),
                "rearranged_electrons_half_l1": density.get("rearranged_electrons_half_l1"),
                "li_density_change_electrons": density.get("li_region", {}).get(
                    "density_change_electrons"
                ),
                "interstitial_accumulation_fraction": density.get(
                    "interstitial_region", {}
                ).get("fraction_of_all_accumulation"),
                "ghost_cavity_accumulation_fraction": density.get(
                    "ghost_cavity_probe", {}
                ).get("fraction_of_all_accumulation"),
            }
        )
    return rows


def run_summary(args: argparse.Namespace) -> int:
    states = [_read_json(path) for path in args.states]
    pairs = [_read_json(path) for path in args.pairs]
    mechanism = _read_json(args.mechanism_summary)
    state_rows = _state_csv_rows(states)
    molecule_rows = _molecule_csv_rows(states)
    pair_rows = _pair_csv_rows(pairs)
    _write_csv(args.states_csv, list(state_rows[0]), state_rows)
    _write_csv(args.molecules_csv, list(molecule_rows[0]), molecule_rows)
    _write_csv(args.pairs_csv, list(pair_rows[0]), pair_rows)
    problems: list[str] = []
    if not mechanism.get("ready"):
        problems.append("Stage-B mechanism smoke summary is not ready")
    if not states or not all(record.get("ready") for record in states):
        problems.append("one or more state-localization analyses are not ready")
    if not pairs or not all(record.get("ready") for record in pairs):
        problems.append("one or more paired density-difference analyses are not ready")
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "campaign": args.campaign,
            "kind": "stage_b_localization",
            "ready": not problems,
            "scientific_status": "NUMERICAL_LOCALIZATION_DIAGNOSTIC_ONLY",
            "problems": problems,
            "interpretation": (
                "This gate validates cube integrity, spin integration, geometric partitioning, "
                "and same-grid density redistribution. Proxy flags require later Multiwfn or "
                "Hirshfeld/Bader cross-checks before chemical interpretation."
            ),
            "source_mechanism_summary": {
                "path": str(Path(args.mechanism_summary).resolve()),
                "sha256": sha256_file(args.mechanism_summary),
            },
            "state_count": len(states),
            "pair_count": len(pairs),
            "states": states,
            "pairs": pairs,
            "tables": {
                "states": str(Path(args.states_csv).resolve()),
                "molecules": str(Path(args.molecules_csv).resolve()),
                "pairs": str(Path(args.pairs_csv).resolve()),
            },
        },
    )
    return 0


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--system", required=True)
    parser.add_argument("--replica", required=True, type=int)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--candidate-metadata", required=True)
    parser.add_argument("--coordinates", required=True)
    parser.add_argument("--cell", required=True)
    parser.add_argument("--spec", required=True)
    parser.add_argument("--methods", required=True)
    parser.add_argument("--systems", required=True)
    parser.add_argument("--output", required=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    state = subparsers.add_parser("state")
    _common(state)
    state.add_argument("--state", required=True)
    state.add_argument("--spin-cube", required=True)
    state.add_argument("--electron-cube", required=True)
    state.set_defaults(func=run_state)
    pair = subparsers.add_parser("pair")
    _common(pair)
    pair.add_argument("--state-records", nargs=2, required=True)
    pair.add_argument("--li0-electron-cube", required=True)
    pair.add_argument("--li-plus-e-electron-cube", required=True)
    pair.add_argument("--difference-cube", required=True)
    pair.set_defaults(func=run_pair)
    summary = subparsers.add_parser("summary")
    summary.add_argument("--campaign", required=True)
    summary.add_argument("--mechanism-summary", required=True)
    summary.add_argument("--states", nargs="+", required=True)
    summary.add_argument("--pairs", nargs="+", required=True)
    summary.add_argument("--states-csv", required=True)
    summary.add_argument("--molecules-csv", required=True)
    summary.add_argument("--pairs-csv", required=True)
    summary.add_argument("--output", required=True)
    summary.set_defaults(func=run_summary)
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
