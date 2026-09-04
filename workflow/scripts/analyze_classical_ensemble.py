#!/usr/bin/env python3
"""Analyze classical trajectories and export deterministic electronic-structure seeds."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from solvelec.cavity_visualization import (
    cavity_centered_positions,
    local_residue_mask,
    render_cavity_hbond_svg,
    write_local_pdb,
    write_pymol_script,
)
from solvelec.composition import AVOGADRO_MOL_INV
from solvelec.hydrogen_bonds import find_hydrogen_bonds, hydrogen_bond_counts
from solvelec.provenance import sha256_file
from solvelec.trajectory import (
    VDW_RADII_ANGSTROM,
    autocorrelation_summary,
    cell_matrix,
    infer_element,
    largest_void_proxy,
    pair_distances,
    select_representative_indices,
)

TIMESERIES_FIELDS = (
    "frame_index",
    "time_ps",
    "elapsed_ps",
    "volume_nm3",
    "density_g_ml",
    "void_radius_angstrom",
    "void_x_angstrom",
    "void_y_angstrom",
    "void_z_angstrom",
    "eda_thf_contacts",
    "eda_thf_hydrogen_bonds",
    "eda_eda_hydrogen_bonds",
    "cavity_associated_hydrogen_bonds",
    "cavity_bridging_hydrogen_bonds",
)

HYDROGEN_BOND_FIELDS = (
    "frame_index",
    "time_ps",
    "elapsed_ps",
    "cavity_x_angstrom",
    "cavity_y_angstrom",
    "cavity_z_angstrom",
    "cavity_radius_angstrom",
    "donor_index",
    "hydrogen_index",
    "acceptor_index",
    "bond_type",
    "donor_acceptor_distance_angstrom",
    "donor_hydrogen_acceptor_angle_degree",
    "cavity_segment_distance_angstrom",
    "cavity_associated",
    "cavity_bridging",
)


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _write_csv(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _element_for_atom(atom: Any) -> str:
    try:
        explicit = str(atom.element)
    except (AttributeError, ValueError):
        explicit = None
    return infer_element(str(atom.name), explicit)


def _atom_groups(universe: Any) -> tuple[list[str], dict[str, np.ndarray]]:
    elements = [_element_for_atom(atom) for atom in universe.atoms]
    groups: dict[str, list[int]] = {
        "thf_o": [],
        "eda_n": [],
        "heavy": [],
    }
    for atom, element in zip(universe.atoms, elements, strict=True):
        residue = str(atom.resname).upper()
        if element != "H":
            groups["heavy"].append(int(atom.index))
        if residue == "THF" and element == "O":
            groups["thf_o"].append(int(atom.index))
        if residue == "EDA" and element == "N":
            groups["eda_n"].append(int(atom.index))
    return elements, {name: np.asarray(indices, dtype=int) for name, indices in groups.items()}


def _bonded_hydrogens(
    universe: Any, nitrogen_indices: np.ndarray, elements: list[str]
) -> dict[int, list[int]]:
    mapping: dict[int, list[int]] = {}
    for index in nitrogen_indices:
        atom = universe.atoms[int(index)]
        try:
            bonded = atom.bonded_atoms
        except (AttributeError, ValueError):
            mapping[int(index)] = []
            continue
        mapping[int(index)] = [
            int(other.index) for other in bonded if elements[int(other.index)] == "H"
        ]
    return mapping


def _rdf_pairs(
    groups: dict[str, np.ndarray],
) -> list[tuple[str, np.ndarray, np.ndarray, bool, bool]]:
    pairs: list[tuple[str, np.ndarray, np.ndarray, bool, bool]] = []
    if len(groups["thf_o"]):
        pairs.append(("thf_o-thf_o", groups["thf_o"], groups["thf_o"], True, False))
    if len(groups["eda_n"]):
        if len(groups["thf_o"]):
            pairs.append(("eda_n-thf_o", groups["eda_n"], groups["thf_o"], False, False))
        pairs.append(("eda_n-eda_n", groups["eda_n"], groups["eda_n"], True, True))
    return pairs


def _trajectory_frame_indices(trajectory: Any, stride_ps: float) -> list[int]:
    if len(trajectory) < 2:
        raise ValueError("trajectory must contain at least two frames")
    timestep_ps = float(trajectory.dt)
    if timestep_ps <= 0:
        raise ValueError("trajectory time step must be positive")
    stride = max(1, int(round(stride_ps / timestep_ps)))
    indices = list(range(0, len(trajectory), stride))
    if indices[-1] != len(trajectory) - 1:
        indices.append(len(trajectory) - 1)
    return indices


def _analyze_trajectory(
    tpr: Path,
    trajectory_path: Path,
    settings: dict[str, Any],
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    import MDAnalysis as mda

    universe = mda.Universe(str(tpr), str(trajectory_path))
    elements, groups = _atom_groups(universe)
    if len(groups["thf_o"]) == 0 and len(groups["eda_n"]) == 0:
        raise ValueError("topology contains neither THF oxygen nor EDA nitrogen atoms")
    heavy_indices = groups["heavy"]
    heavy_radii = np.asarray([VDW_RADII_ANGSTROM[elements[index]] for index in heavy_indices])
    nitrogen_to_hydrogen = _bonded_hydrogens(universe, groups["eda_n"], elements)
    residue_ids = np.asarray(universe.atoms.resindices, dtype=int)
    frame_indices = _trajectory_frame_indices(
        universe.trajectory, float(settings["analysis_stride_ps"])
    )
    rdf_max = float(settings["rdf_max_angstrom"])
    rdf_width = float(settings["rdf_bin_width_angstrom"])
    rdf_edges = np.arange(0.0, rdf_max + rdf_width * 0.5, rdf_width)
    if len(rdf_edges) < 2 or rdf_edges[-1] < rdf_max * 0.99:
        raise ValueError("RDF bin settings do not cover the requested range")
    rdf_counts = {
        name: np.zeros(len(rdf_edges) - 1, dtype=np.int64)
        for name, _first, _second, _same, _exclude in _rdf_pairs(groups)
    }
    rdf_expected = {
        name: np.zeros(len(rdf_edges) - 1, dtype=float)
        for name, _first, _second, _same, _exclude in _rdf_pairs(groups)
    }
    shell_volumes = (4.0 / 3.0) * np.pi * (rdf_edges[1:] ** 3 - rdf_edges[:-1] ** 3)
    total_mass_g_mol = float(universe.atoms.total_mass())
    first_time_ps = float(universe.trajectory[frame_indices[0]].time)
    rows: list[dict[str, Any]] = []
    hydrogen_bond_rows: list[dict[str, Any]] = []
    for frame_index in frame_indices:
        frame = universe.trajectory[frame_index]
        matrix = cell_matrix(frame.dimensions)
        shortest_vector = min(float(np.linalg.norm(vector)) for vector in matrix)
        if rdf_max >= 0.5 * shortest_vector:
            raise ValueError(f"RDF maximum {rdf_max:g} A is not below half the shortest box vector")
        positions = np.asarray(universe.atoms.positions, dtype=float)
        volume_angstrom3 = abs(float(np.linalg.det(matrix)))
        volume_nm3 = volume_angstrom3 / 1000.0
        density = total_mass_g_mol / AVOGADRO_MOL_INV / (volume_nm3 * 1.0e-21)
        void = largest_void_proxy(
            positions[heavy_indices],
            heavy_radii,
            matrix,
            points_per_axis=int(settings["void_grid_points_per_axis"]),
            refinement_levels=int(settings["void_refinement_levels"]),
        )
        contacts = 0
        bond_counts = hydrogen_bond_counts([])
        if len(groups["eda_n"]):
            if len(groups["thf_o"]):
                contact_distances = pair_distances(
                    positions[groups["eda_n"]], positions[groups["thf_o"]], matrix
                )
                contacts = int(
                    np.count_nonzero(
                        contact_distances <= float(settings["eda_thf_contact_cutoff_angstrom"])
                    )
                )
            bonds = find_hydrogen_bonds(
                positions,
                nitrogen_to_hydrogen,
                {"thf_o": groups["thf_o"], "eda_n": groups["eda_n"]},
                residue_ids,
                matrix,
                np.asarray(void["cartesian_angstrom"], dtype=float),
                float(void["radius_angstrom"]),
                distance_cutoff_angstrom=float(settings["hydrogen_bond_distance_angstrom"]),
                angle_cutoff_degree=float(settings["hydrogen_bond_angle_degree"]),
                cavity_shell_thickness_angstrom=float(
                    settings["cavity_hbond_shell_thickness_angstrom"]
                ),
                cavity_bridge_margin_angstrom=float(
                    settings["cavity_hbond_bridge_margin_angstrom"]
                ),
            )
            bond_counts = hydrogen_bond_counts(bonds)
            for bond in bonds:
                hydrogen_bond_rows.append(
                    {
                        "frame_index": int(frame_index),
                        "time_ps": float(frame.time),
                        "elapsed_ps": float(frame.time) - first_time_ps,
                        "cavity_x_angstrom": float(void["cartesian_angstrom"][0]),
                        "cavity_y_angstrom": float(void["cartesian_angstrom"][1]),
                        "cavity_z_angstrom": float(void["cartesian_angstrom"][2]),
                        "cavity_radius_angstrom": float(void["radius_angstrom"]),
                        **bond.as_dict(),
                    }
                )
        center = void["cartesian_angstrom"]
        rows.append(
            {
                "frame_index": int(frame_index),
                "time_ps": float(frame.time),
                "elapsed_ps": float(frame.time) - first_time_ps,
                "volume_nm3": volume_nm3,
                "density_g_ml": density,
                "void_radius_angstrom": float(void["radius_angstrom"]),
                "void_x_angstrom": float(center[0]),
                "void_y_angstrom": float(center[1]),
                "void_z_angstrom": float(center[2]),
                "eda_thf_contacts": contacts,
                **bond_counts,
            }
        )
        for name, first_indices, second_indices, same, exclude in _rdf_pairs(groups):
            distances = pair_distances(
                positions[first_indices],
                positions[second_indices],
                matrix,
                same_group=same,
                first_residue_ids=np.asarray(universe.atoms[first_indices].resindices),
                second_residue_ids=np.asarray(universe.atoms[second_indices].resindices),
                exclude_same_residue=exclude,
            )
            rdf_counts[name] += np.histogram(distances, bins=rdf_edges)[0]
            rdf_expected[name] += len(distances) * shell_volumes / volume_angstrom3

    rdf_rows: list[dict[str, Any]] = []
    for name in rdf_counts:
        values = np.divide(
            rdf_counts[name],
            rdf_expected[name],
            out=np.zeros_like(rdf_expected[name]),
            where=rdf_expected[name] > 0,
        )
        for index, value in enumerate(values):
            rdf_rows.append(
                {
                    "pair": name,
                    "r_lower_angstrom": float(rdf_edges[index]),
                    "r_upper_angstrom": float(rdf_edges[index + 1]),
                    "r_center_angstrom": float((rdf_edges[index] + rdf_edges[index + 1]) / 2),
                    "g_r": float(value),
                    "count": int(rdf_counts[name][index]),
                    "expected_count": float(rdf_expected[name][index]),
                }
            )

    autocorrelation: dict[str, Any] = {}
    times = [float(row["elapsed_ps"]) for row in rows]
    descriptors = ["volume_nm3", "density_g_ml", "void_radius_angstrom"]
    if len(groups["eda_n"]):
        descriptors.extend(
            [
                "eda_eda_hydrogen_bonds",
                "cavity_associated_hydrogen_bonds",
                "cavity_bridging_hydrogen_bonds",
            ]
        )
        if len(groups["thf_o"]):
            descriptors.extend(["eda_thf_contacts", "eda_thf_hydrogen_bonds"])
    for descriptor in descriptors:
        autocorrelation[descriptor] = autocorrelation_summary(
            times, [float(row[descriptor]) for row in rows]
        )
    minimum_effective_samples = float(settings["minimum_effective_samples"])
    checks = {
        "minimum_analysis_frames": len(rows) >= int(settings["minimum_analysis_frames"]),
        "effective_samples": all(
            float(metrics["effective_sample_size"]) >= minimum_effective_samples
            for metrics in autocorrelation.values()
        ),
        "positive_void_proxy": all(float(row["void_radius_angstrom"]) > 0 for row in rows),
    }
    metadata = {
        "ready": all(checks.values()),
        "checks": checks,
        "frame_count": len(rows),
        "trajectory_frame_count": len(universe.trajectory),
        "analysis_stride_ps": float(settings["analysis_stride_ps"]),
        "sampled_duration_ns": (float(rows[-1]["elapsed_ps"]) - float(rows[0]["elapsed_ps"]))
        / 1000.0,
        "atom_count": len(universe.atoms),
        "atom_groups": {name: int(len(indices)) for name, indices in groups.items()},
        "hydrogen_bond_donor_count": sum(bool(value) for value in nitrogen_to_hydrogen.values()),
        "hydrogen_bond_definition": {
            "distance_cutoff_angstrom": float(settings["hydrogen_bond_distance_angstrom"]),
            "angle_cutoff_degree": float(settings["hydrogen_bond_angle_degree"]),
            "cavity_shell_thickness_angstrom": float(
                settings["cavity_hbond_shell_thickness_angstrom"]
            ),
            "cavity_bridge_margin_angstrom": float(settings["cavity_hbond_bridge_margin_angstrom"]),
            "interpretation": "geometric cavity relation; no bond to the cavity center is implied",
        },
        "autocorrelation": autocorrelation,
        "mean_descriptors": {
            descriptor: float(np.mean([float(row[descriptor]) for row in rows]))
            for descriptor in descriptors
        },
    }
    return metadata, rows, rdf_rows, hydrogen_bond_rows


def run_analyze(args: argparse.Namespace) -> int:
    spec_path = Path(args.spec)
    methods_path = Path(args.methods)
    validation_path = Path(args.classical_validation)
    tpr_path = Path(args.tpr)
    trajectory_path = Path(args.trajectory)
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    methods = json.loads(methods_path.read_text(encoding="utf-8"))
    classical_validation = json.loads(validation_path.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    rdf_rows: list[dict[str, Any]] = []
    hydrogen_bond_rows: list[dict[str, Any]] = []
    try:
        metrics, rows, rdf_rows, hydrogen_bond_rows = _analyze_trajectory(
            tpr_path, trajectory_path, methods["trajectory_analysis"]
        )
        metrics["checks"]["classical_pilot_ready"] = bool(
            classical_validation.get("metrics", {}).get("ready")
        )
        metrics["ready"] = all(metrics["checks"].values())
        result = {
            "schema_version": 1,
            "system_id": spec["system_id"],
            "replica": int(spec["replica"]),
            "ready": bool(metrics["ready"]),
            "metrics": metrics,
            "inputs": {
                "spec": {"path": str(spec_path.resolve()), "sha256": sha256_file(spec_path)},
                "tpr": {"path": str(tpr_path.resolve()), "sha256": sha256_file(tpr_path)},
                "trajectory": {
                    "path": str(trajectory_path.resolve()),
                    "sha256": sha256_file(trajectory_path),
                },
                "classical_validation": {
                    "path": str(validation_path.resolve()),
                    "sha256": sha256_file(validation_path),
                },
            },
        }
    except Exception as exc:
        result = {
            "schema_version": 1,
            "system_id": spec.get("system_id"),
            "replica": spec.get("replica"),
            "ready": False,
            "metrics": {"ready": False, "error": str(exc)},
        }
    _write_csv(Path(args.timeseries), TIMESERIES_FIELDS, rows)
    _write_csv(
        Path(args.rdf),
        (
            "pair",
            "r_lower_angstrom",
            "r_upper_angstrom",
            "r_center_angstrom",
            "g_r",
            "count",
            "expected_count",
        ),
        rdf_rows,
    )
    _write_csv(Path(args.hydrogen_bonds), HYDROGEN_BOND_FIELDS, hydrogen_bond_rows)
    result["outputs"] = {
        "timeseries": str(Path(args.timeseries).resolve()),
        "rdf": str(Path(args.rdf).resolve()),
        "hydrogen_bonds": str(Path(args.hydrogen_bonds).resolve()),
    }
    _write_json(Path(args.output), result)
    return 0


def _read_timeseries(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append({name: float(value) for name, value in row.items() if value is not None})
    return rows


def _write_snapshot(
    universe: Any,
    frame_index: int,
    xyz_path: Path,
    cell_path: Path,
    comment: str,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    frame = universe.trajectory[frame_index]
    matrix = cell_matrix(frame.dimensions)
    elements = [_element_for_atom(atom) for atom in universe.atoms]
    try:
        positions = np.asarray(
            universe.atoms.wrap(compound="residues", center="com", inplace=False), dtype=float
        )
    except (AttributeError, ValueError):
        positions = np.asarray(universe.atoms.wrap(compound="atoms", inplace=False), dtype=float)
    lines = [str(len(elements)), comment]
    for element, position in zip(elements, positions, strict=True):
        label = "Li" if element == "LI" else element
        lines.append(f"{label:<2s} {position[0]: .10f} {position[1]: .10f} {position[2]: .10f}")
    xyz_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_xyz = xyz_path.with_name(f".{xyz_path.name}.tmp")
    temporary_xyz.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary_xyz.replace(xyz_path)
    cell_lines = ["&CELL"]
    for name, vector in zip(("A", "B", "C"), matrix, strict=True):
        cell_lines.append(f"  {name} {vector[0]:.10f} {vector[1]:.10f} {vector[2]:.10f}")
    cell_lines.extend(["  PERIODIC XYZ", "&END CELL"])
    cell_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_cell = cell_path.with_name(f".{cell_path.name}.tmp")
    temporary_cell.write_text("\n".join(cell_lines) + "\n", encoding="utf-8")
    temporary_cell.replace(cell_path)
    return elements, positions, matrix


def run_select(args: argparse.Namespace) -> int:
    import MDAnalysis as mda

    analysis_path = Path(args.analysis)
    methods_path = Path(args.methods)
    tpr_path = Path(args.tpr)
    trajectory_path = Path(args.trajectory)
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    methods = json.loads(methods_path.read_text(encoding="utf-8"))
    settings = methods["trajectory_analysis"]
    if not analysis.get("ready"):
        raise ValueError(f"trajectory analysis is not ready: {analysis_path}")
    snapshot_count = int(settings["snapshots_per_replica"])
    if snapshot_count != 1:
        raise ValueError("stage-A snapshot export currently requires snapshots_per_replica = 1")
    rows = _read_timeseries(Path(args.timeseries))
    feature_names = ["density_g_ml", "void_radius_angstrom"]
    for candidate in (
        "eda_thf_contacts",
        "eda_thf_hydrogen_bonds",
        "eda_eda_hydrogen_bonds",
        "cavity_associated_hydrogen_bonds",
        "cavity_bridging_hydrogen_bonds",
    ):
        if rows and candidate in rows[0] and candidate in analysis["metrics"]["autocorrelation"]:
            feature_names.append(candidate)
    max_tau_ps = max(
        float(record["integrated_autocorrelation_time_ps"])
        for record in analysis["metrics"]["autocorrelation"].values()
    )
    separation_ps = max(
        float(settings["minimum_snapshot_separation_ps"]),
        float(settings["decorrelation_multiplier"]) * max_tau_ps,
    )
    selection_rows = [{**row, "time_ps": float(row["elapsed_ps"])} for row in rows]
    selected_index = select_representative_indices(
        selection_rows,
        feature_names,
        count=1,
        minimum_time_ps=float(settings["equilibrated_start_ns"]) * 1000.0,
        minimum_separation_ps=separation_ps,
    )[0]
    selected = rows[selected_index]
    universe = mda.Universe(str(tpr_path), str(trajectory_path))
    snapshot_id = (
        f"{analysis['system_id']}_r{int(analysis['replica'])}_"
        f"t{float(selected['elapsed_ps']):.0f}ps"
    )
    elements, positions, matrix = _write_snapshot(
        universe,
        int(selected["frame_index"]),
        Path(args.xyz),
        Path(args.cell),
        f"snapshot_id={snapshot_id} source_time_ps={float(selected['time_ps']):.6f}",
    )
    _, groups = _atom_groups(universe)
    residue_ids = np.asarray(universe.atoms.resindices, dtype=int)
    cavity_center = np.asarray(
        [
            selected["void_x_angstrom"],
            selected["void_y_angstrom"],
            selected["void_z_angstrom"],
        ],
        dtype=float,
    )
    cavity_radius = float(selected["void_radius_angstrom"])
    donor_to_hydrogens = _bonded_hydrogens(universe, groups["eda_n"], elements)
    snapshot_bonds = find_hydrogen_bonds(
        positions,
        donor_to_hydrogens,
        {"thf_o": groups["thf_o"], "eda_n": groups["eda_n"]},
        residue_ids,
        matrix,
        cavity_center,
        cavity_radius,
        distance_cutoff_angstrom=float(settings["hydrogen_bond_distance_angstrom"]),
        angle_cutoff_degree=float(settings["hydrogen_bond_angle_degree"]),
        cavity_shell_thickness_angstrom=float(settings["cavity_hbond_shell_thickness_angstrom"]),
        cavity_bridge_margin_angstrom=float(settings["cavity_hbond_bridge_margin_angstrom"]),
    )
    centered_positions = cavity_centered_positions(
        positions, residue_ids, groups["heavy"], cavity_center, matrix
    )
    visualization_radius = float(settings["cavity_visualization_radius_angstrom"])
    local_mask = local_residue_mask(
        centered_positions, residue_ids, groups["heavy"], visualization_radius
    )
    local_indices = set(int(value) for value in np.flatnonzero(local_mask))
    local_atoms = [
        {
            "index": int(atom.index),
            "name": str(atom.name),
            "resname": str(atom.resname),
            "resid": int(atom.resid),
            "element": elements[int(atom.index)],
            "position_angstrom": centered_positions[int(atom.index)].tolist(),
        }
        for atom in universe.atoms
        if int(atom.index) in local_indices
    ]
    try:
        covalent_bonds = [
            (int(bond.atoms[0].index), int(bond.atoms[1].index))
            for bond in universe.bonds
            if int(bond.atoms[0].index) in local_indices
            and int(bond.atoms[1].index) in local_indices
        ]
    except (AttributeError, ValueError):
        covalent_bonds = []
    local_bonds = [
        bond
        for bond in snapshot_bonds
        if bond.donor_index in local_indices
        and bond.hydrogen_index in local_indices
        and bond.acceptor_index in local_indices
    ]
    local_hydrogen_bonds = [bond.as_dict() for bond in local_bonds]
    local_pdb_path = Path(args.local_pdb)
    serials = write_local_pdb(local_pdb_path, local_atoms)
    write_pymol_script(
        Path(args.pymol),
        local_pdb_path,
        cavity_radius,
        local_hydrogen_bonds,
        serials,
    )
    render_cavity_hbond_svg(
        Path(args.svg),
        title=f"{snapshot_id}: cavity-centered hydrogen-bond shell",
        atoms=local_atoms,
        covalent_bonds=covalent_bonds,
        hydrogen_bonds=local_hydrogen_bonds,
        cavity_radius_angstrom=cavity_radius,
        view_radius_angstrom=visualization_radius,
    )
    cavity_hbond_record = {
        "schema_version": 1,
        "system_id": analysis["system_id"],
        "replica": int(analysis["replica"]),
        "snapshot_id": snapshot_id,
        "scientific_status": "GEOMETRIC_CAVITY_HBOND_VISUALIZATION",
        "interpretation": (
            "Hydrogen bonds are classified relative to a geometric void proxy; "
            "no bond to the cavity center or electron localization is implied."
        ),
        "cavity": {
            "center_in_periodic_cell_angstrom": cavity_center.tolist(),
            "center_in_local_view_angstrom": [0.0, 0.0, 0.0],
            "radius_angstrom": cavity_radius,
            "visualization_radius_angstrom": visualization_radius,
        },
        "hydrogen_bond_definition": analysis["metrics"]["hydrogen_bond_definition"],
        "hydrogen_bond_counts_whole_snapshot": hydrogen_bond_counts(snapshot_bonds),
        "hydrogen_bond_counts_local_view": hydrogen_bond_counts(local_bonds),
        "local_atom_count": len(local_atoms),
        "local_hydrogen_bonds": local_hydrogen_bonds,
        "outputs": {
            "local_pdb": str(local_pdb_path.resolve()),
            "pymol": str(Path(args.pymol).resolve()),
            "svg": str(Path(args.svg).resolve()),
        },
    }
    _write_json(Path(args.cavity_hbonds), cavity_hbond_record)
    metadata = {
        "schema_version": 1,
        "ready": True,
        "system_id": analysis["system_id"],
        "replica": int(analysis["replica"]),
        "snapshot_id": snapshot_id,
        "li_atom_present": False,
        "selection": {
            "method": "equilibrated robust-medoid",
            "features": feature_names,
            "candidate_start_ns": float(settings["equilibrated_start_ns"]),
            "required_separation_ps": separation_ps,
            "frame_index": int(selected["frame_index"]),
            "source_time_ps": float(selected["time_ps"]),
            "elapsed_ps": float(selected["elapsed_ps"]),
            "descriptors": {name: float(selected[name]) for name in feature_names},
        },
        "structure": {
            "atom_count": len(elements),
            "cell_vectors_angstrom": matrix.tolist(),
            "xyz": {"path": str(Path(args.xyz).resolve()), "sha256": sha256_file(args.xyz)},
            "cell": {"path": str(Path(args.cell).resolve()), "sha256": sha256_file(args.cell)},
            "cavity_hbonds": {
                "path": str(Path(args.cavity_hbonds).resolve()),
                "sha256": sha256_file(args.cavity_hbonds),
            },
            "cavity_svg": {
                "path": str(Path(args.svg).resolve()),
                "sha256": sha256_file(args.svg),
            },
            "cavity_local_pdb": {
                "path": str(local_pdb_path.resolve()),
                "sha256": sha256_file(local_pdb_path),
            },
            "cavity_pymol": {
                "path": str(Path(args.pymol).resolve()),
                "sha256": sha256_file(args.pymol),
            },
        },
        "source": {
            "analysis": {
                "path": str(analysis_path.resolve()),
                "sha256": sha256_file(analysis_path),
            },
            "tpr_sha256": analysis["inputs"]["tpr"]["sha256"],
            "trajectory_sha256": analysis["inputs"]["trajectory"]["sha256"],
        },
    }
    _write_json(Path(args.output), metadata)
    return 0


def summarize_records(records: list[dict[str, Any]], kind: str) -> dict[str, Any]:
    if not records:
        raise ValueError("at least one record is required")
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[str(record["system_id"])].append(record)
    systems: dict[str, Any] = {}
    for system_id, group in sorted(grouped.items()):
        replicas = [int(record["replica"]) for record in group]
        checks = {
            "all_records_ready": all(bool(record.get("ready")) for record in group),
            "unique_replicas": len(replicas) == len(set(replicas)),
        }
        system_result: dict[str, Any] = {
            "ready": all(checks.values()),
            "checks": checks,
            "replicas": sorted(replicas),
        }
        if kind == "analysis":
            descriptor_maps = [
                record.get("metrics", {}).get("mean_descriptors", {}) for record in group
            ]
            common_descriptors = (
                set.intersection(*(set(values) for values in descriptor_maps))
                if descriptor_maps
                else set()
            )
            ensemble_descriptors: dict[str, Any] = {}
            for descriptor in sorted(common_descriptors):
                replica_means = [float(values[descriptor]) for values in descriptor_maps]
                ensemble_descriptors[descriptor] = {
                    "mean": float(np.mean(replica_means)),
                    "replica_standard_deviation": (
                        float(np.std(replica_means, ddof=1)) if len(replica_means) > 1 else 0.0
                    ),
                    "replica_means": replica_means,
                }
            effective_samples = [
                float(value["effective_sample_size"])
                for record in group
                for value in record.get("metrics", {}).get("autocorrelation", {}).values()
                if "effective_sample_size" in value
            ]
            void_radii = [
                float(record["metrics"]["mean_descriptors"]["void_radius_angstrom"])
                for record in group
                if "mean_descriptors" in record.get("metrics", {})
                and "void_radius_angstrom" in record["metrics"]["mean_descriptors"]
            ]
            system_result["minimum_effective_samples"] = (
                min(effective_samples) if effective_samples else None
            )
            system_result["mean_void_radius_angstrom"] = (
                float(np.mean(void_radii)) if len(void_radii) == len(group) else None
            )
            system_result["ensemble_mean_descriptors"] = ensemble_descriptors
        systems[system_id] = system_result
    return {
        "schema_version": 1,
        "kind": kind,
        "ready": all(result["ready"] for result in systems.values()),
        "record_count": len(records),
        "systems": systems,
    }


def run_summary(args: argparse.Namespace) -> int:
    records = [json.loads(Path(path).read_text(encoding="utf-8")) for path in args.inputs]
    result = summarize_records(records, args.kind)
    result["campaign"] = args.campaign
    _write_json(Path(args.output), result)
    return 0


def run_gate(args: argparse.Namespace) -> int:
    summary_path = Path(args.summary)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    output = Path(args.output)
    output.unlink(missing_ok=True)
    if not summary.get("ready"):
        print(f"ERROR: stage-A analysis failed; inspect {summary_path}", file=sys.stderr)
        return 4
    digest = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(f"summary_sha256={digest}\n", encoding="utf-8")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser("analyze")
    analyze.add_argument("--spec", required=True)
    analyze.add_argument("--methods", required=True)
    analyze.add_argument("--classical-validation", required=True)
    analyze.add_argument("--tpr", required=True)
    analyze.add_argument("--trajectory", required=True)
    analyze.add_argument("--timeseries", required=True)
    analyze.add_argument("--rdf", required=True)
    analyze.add_argument("--hydrogen-bonds", required=True)
    analyze.add_argument("--output", required=True)
    analyze.set_defaults(func=run_analyze)

    select = subparsers.add_parser("select")
    select.add_argument("--analysis", required=True)
    select.add_argument("--timeseries", required=True)
    select.add_argument("--methods", required=True)
    select.add_argument("--tpr", required=True)
    select.add_argument("--trajectory", required=True)
    select.add_argument("--xyz", required=True)
    select.add_argument("--cell", required=True)
    select.add_argument("--cavity-hbonds", required=True)
    select.add_argument("--local-pdb", required=True)
    select.add_argument("--pymol", required=True)
    select.add_argument("--svg", required=True)
    select.add_argument("--output", required=True)
    select.set_defaults(func=run_select)

    summary = subparsers.add_parser("summary")
    summary.add_argument("--campaign", required=True)
    summary.add_argument("--kind", choices=("analysis", "snapshot"), required=True)
    summary.add_argument("--output", required=True)
    summary.add_argument("inputs", nargs="+")
    summary.set_defaults(func=run_summary)

    gate = subparsers.add_parser("gate")
    gate.add_argument("--summary", required=True)
    gate.add_argument("--output", required=True)
    gate.set_defaults(func=run_gate)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
