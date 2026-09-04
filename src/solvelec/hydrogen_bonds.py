from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .trajectory import minimum_image_vectors


@dataclass(frozen=True)
class HydrogenBond:
    """One geometrically identified N-H...acceptor interaction."""

    donor_index: int
    hydrogen_index: int
    acceptor_index: int
    bond_type: str
    donor_acceptor_distance_angstrom: float
    donor_hydrogen_acceptor_angle_degree: float
    cavity_segment_distance_angstrom: float
    cavity_associated: bool
    cavity_bridging: bool

    def as_dict(self) -> dict[str, int | float | str | bool]:
        return asdict(self)


def point_segment_distance(point: np.ndarray, start: np.ndarray, end: np.ndarray) -> float:
    """Return the shortest Euclidean distance from a point to a finite segment."""

    point_array = np.asarray(point, dtype=float)
    start_array = np.asarray(start, dtype=float)
    end_array = np.asarray(end, dtype=float)
    if point_array.shape != (3,) or start_array.shape != (3,) or end_array.shape != (3,):
        raise ValueError("point and segment endpoints must be three-vectors")
    segment = end_array - start_array
    squared_length = float(segment @ segment)
    if squared_length == 0:
        return float(np.linalg.norm(point_array - start_array))
    fraction = float((point_array - start_array) @ segment / squared_length)
    fraction = min(1.0, max(0.0, fraction))
    closest = start_array + fraction * segment
    return float(np.linalg.norm(point_array - closest))


def find_hydrogen_bonds(
    positions: np.ndarray,
    donor_to_hydrogens: dict[int, list[int]],
    acceptor_groups: dict[str, np.ndarray],
    residue_ids: np.ndarray,
    cell: np.ndarray,
    cavity_center: np.ndarray,
    cavity_radius_angstrom: float,
    *,
    distance_cutoff_angstrom: float,
    angle_cutoff_degree: float,
    cavity_shell_thickness_angstrom: float,
    cavity_bridge_margin_angstrom: float,
) -> list[HydrogenBond]:
    """Find EDA hydrogen bonds and classify their relation to a geometric void.

    ``cavity_associated`` means that the donor--acceptor segment approaches the
    void surface within the configured shell thickness. ``cavity_bridging`` is
    the stricter geometric test that the segment intersects the void sphere,
    allowing only a small numerical margin. Neither label implies bonding to an
    electron or to the cavity center.
    """

    coordinates = np.asarray(positions, dtype=float)
    residues = np.asarray(residue_ids)
    matrix = np.asarray(cell, dtype=float)
    center = np.asarray(cavity_center, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("positions must be shaped (n, 3)")
    if residues.shape != (len(coordinates),):
        raise ValueError("one residue id is required per atom")
    if center.shape != (3,) or cavity_radius_angstrom <= 0:
        raise ValueError("a three-vector cavity center and positive radius are required")
    if distance_cutoff_angstrom <= 0 or not 0 < angle_cutoff_degree <= 180:
        raise ValueError("hydrogen-bond distance and angle thresholds are invalid")
    if cavity_shell_thickness_angstrom < 0 or cavity_bridge_margin_angstrom < 0:
        raise ValueError("cavity shell and bridge margins must be non-negative")

    bonds: list[HydrogenBond] = []
    for donor_index in sorted(donor_to_hydrogens):
        hydrogen_indices = sorted(donor_to_hydrogens[donor_index])
        if not hydrogen_indices:
            continue
        donor_position = coordinates[[donor_index]]
        donor_local = minimum_image_vectors(center[None, :], donor_position, matrix)[0, 0]
        for acceptor_kind, raw_indices in sorted(acceptor_groups.items()):
            acceptor_indices = np.asarray(raw_indices, dtype=int)
            if not len(acceptor_indices):
                continue
            allowed = acceptor_indices != donor_index
            if acceptor_kind == "eda_n":
                allowed &= residues[acceptor_indices] != residues[donor_index]
            acceptor_indices = acceptor_indices[allowed]
            if not len(acceptor_indices):
                continue
            acceptor_positions = coordinates[acceptor_indices]
            donor_to_acceptor = minimum_image_vectors(donor_position, acceptor_positions, matrix)[0]
            distances = np.linalg.norm(donor_to_acceptor, axis=1)
            nearby_offsets = np.flatnonzero(distances <= distance_cutoff_angstrom)
            for hydrogen_index in hydrogen_indices:
                hydrogen_position = coordinates[[hydrogen_index]]
                hydrogen_to_donor = minimum_image_vectors(
                    hydrogen_position, donor_position, matrix
                )[0, 0]
                hydrogen_to_acceptor = minimum_image_vectors(
                    hydrogen_position,
                    acceptor_positions[nearby_offsets],
                    matrix,
                )[0]
                denominator = np.linalg.norm(hydrogen_to_donor) * np.linalg.norm(
                    hydrogen_to_acceptor, axis=1
                )
                valid = denominator > 0
                cosines = np.ones(len(nearby_offsets), dtype=float)
                cosines[valid] = (
                    hydrogen_to_acceptor[valid] @ hydrogen_to_donor / denominator[valid]
                )
                angles = np.rad2deg(np.arccos(np.clip(cosines, -1.0, 1.0)))
                for local_offset in np.flatnonzero(angles >= angle_cutoff_degree):
                    nearby_offset = int(nearby_offsets[local_offset])
                    acceptor_index = int(acceptor_indices[nearby_offset])
                    acceptor_local = donor_local + donor_to_acceptor[nearby_offset]
                    segment_distance = point_segment_distance(
                        np.zeros(3), donor_local, acceptor_local
                    )
                    bonds.append(
                        HydrogenBond(
                            donor_index=int(donor_index),
                            hydrogen_index=int(hydrogen_index),
                            acceptor_index=acceptor_index,
                            bond_type=f"eda_nh_{acceptor_kind}",
                            donor_acceptor_distance_angstrom=float(distances[nearby_offset]),
                            donor_hydrogen_acceptor_angle_degree=float(angles[local_offset]),
                            cavity_segment_distance_angstrom=segment_distance,
                            cavity_associated=segment_distance
                            <= cavity_radius_angstrom + cavity_shell_thickness_angstrom,
                            cavity_bridging=segment_distance
                            <= cavity_radius_angstrom + cavity_bridge_margin_angstrom,
                        )
                    )
    return sorted(
        bonds,
        key=lambda item: (
            item.donor_index,
            item.hydrogen_index,
            item.acceptor_index,
            item.bond_type,
        ),
    )


def hydrogen_bond_counts(bonds: list[HydrogenBond]) -> dict[str, int]:
    counts = {
        "eda_thf_hydrogen_bonds": 0,
        "eda_eda_hydrogen_bonds": 0,
        "cavity_associated_hydrogen_bonds": 0,
        "cavity_bridging_hydrogen_bonds": 0,
    }
    for bond in bonds:
        if bond.bond_type == "eda_nh_thf_o":
            counts["eda_thf_hydrogen_bonds"] += 1
        elif bond.bond_type == "eda_nh_eda_n":
            counts["eda_eda_hydrogen_bonds"] += 1
        if bond.cavity_associated:
            counts["cavity_associated_hydrogen_bonds"] += 1
        if bond.cavity_bridging:
            counts["cavity_bridging_hydrogen_bonds"] += 1
    return counts
