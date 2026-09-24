from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .localization_analysis import MolecularTopology


def minimum_image_vector(
    point: Sequence[float],
    reference: Sequence[float],
    cell_angstrom: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Return the periodic vector from ``reference`` to ``point``."""

    cell = np.asarray(cell_angstrom, dtype=float)
    inverse = np.linalg.inv(cell)
    fractional = (
        np.asarray(point, dtype=float) - np.asarray(reference, dtype=float)
    ) @ inverse
    fractional -= np.rint(fractional)
    return fractional @ cell


def infer_nh_donors(
    elements: Sequence[str],
    positions_angstrom: NDArray[np.float64],
    cell_angstrom: NDArray[np.float64],
    topology: MolecularTopology,
    *,
    maximum_nh_bond_angstrom: float,
) -> tuple[dict[str, Any], ...]:
    """Infer covalent N-H pairs inside already inferred solvent molecules.

    A hydrogen is assigned to at most one nitrogen: the nearest nitrogen in the
    same molecule, provided that the periodic distance is below the configured
    covalent cutoff.  This distinguishes primary/secondary amines from the
    N-H-free TMEDA control without relying on molecule names.
    """

    if maximum_nh_bond_angstrom <= 0:
        raise ValueError("maximum N-H bond distance must be positive")
    positions = np.asarray(positions_angstrom, dtype=float)
    if positions.shape != (len(elements), 3):
        raise ValueError("element and position arrays are inconsistent")

    normalized = [str(element).upper() for element in elements]
    donors: list[dict[str, Any]] = []
    for molecule in topology.molecules:
        atom_indices = [int(value) - 1 for value in molecule["atom_indices_1based"]]
        nitrogen = [index for index in atom_indices if normalized[index] == "N"]
        hydrogen = [index for index in atom_indices if normalized[index] == "H"]
        for hydrogen_index in hydrogen:
            candidates = [
                (
                    float(
                        np.linalg.norm(
                            minimum_image_vector(
                                positions[hydrogen_index],
                                positions[nitrogen_index],
                                cell_angstrom,
                            )
                        )
                    ),
                    nitrogen_index,
                )
                for nitrogen_index in nitrogen
            ]
            if not candidates:
                continue
            distance, nitrogen_index = min(candidates)
            if distance > maximum_nh_bond_angstrom:
                continue
            donors.append(
                {
                    "donor_id": len(donors) + 1,
                    "molecule_id": int(molecule["molecule_id"]),
                    "component": molecule.get("component"),
                    "nitrogen_atom_index_1based": nitrogen_index + 1,
                    "hydrogen_atom_index_1based": hydrogen_index + 1,
                    "nh_bond_length_angstrom": distance,
                }
            )
    return tuple(donors)


def _angle_degrees(first: NDArray[np.float64], second: NDArray[np.float64]) -> float:
    denominator = float(np.linalg.norm(first) * np.linalg.norm(second))
    if denominator <= 0:
        raise ValueError("cannot evaluate an angle involving a zero-length vector")
    cosine = float(np.dot(first, second) / denominator)
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def nh_environment(
    point_angstrom: Sequence[float],
    donors: Sequence[Mapping[str, Any]],
    positions_angstrom: NDArray[np.float64],
    cell_angstrom: NDArray[np.float64],
    *,
    shell_radii_angstrom: Sequence[float],
    inward_angle_minimum_degrees: float,
    distance_decay_angstrom: float,
) -> dict[str, Any]:
    """Describe N-H groups around a point without claiming a hydrogen bond.

    ``N-H-center`` is the conventional angle at H (H->N versus H->center), so
    180 degrees means that the N-H bond points directly toward the center.
    The smooth score is only a deterministic ranking aid for selecting existing
    solvent fluctuations; it is not an interaction energy.
    """

    radii = sorted({float(value) for value in shell_radii_angstrom})
    if not radii or radii[0] <= 0:
        raise ValueError("N-H shell radii must be positive")
    if not 90.0 <= inward_angle_minimum_degrees <= 180.0:
        raise ValueError("inward N-H angle threshold must be between 90 and 180 degrees")
    if distance_decay_angstrom <= 0:
        raise ValueError("N-H distance decay must be positive")

    point = np.asarray(point_angstrom, dtype=float)
    positions = np.asarray(positions_angstrom, dtype=float)
    records: list[dict[str, Any]] = []
    for donor in donors:
        nitrogen_index = int(donor["nitrogen_atom_index_1based"]) - 1
        hydrogen_index = int(donor["hydrogen_atom_index_1based"]) - 1
        h_to_n = minimum_image_vector(
            positions[nitrogen_index], positions[hydrogen_index], cell_angstrom
        )
        h_to_center = minimum_image_vector(
            point, positions[hydrogen_index], cell_angstrom
        )
        distance = float(np.linalg.norm(h_to_center))
        angle = _angle_degrees(h_to_n, h_to_center)
        inward = angle >= inward_angle_minimum_degrees
        angular_weight = max(0.0, (angle - 90.0) / 90.0)
        distance_weight = float(np.exp(-distance / distance_decay_angstrom))
        records.append(
            {
                **dict(donor),
                "hydrogen_to_center_angstrom": distance,
                "n_h_center_angle_degrees": angle,
                "inward_facing": inward,
                "ranking_weight": distance_weight * angular_weight,
            }
        )
    records.sort(
        key=lambda record: (
            float(record["hydrogen_to_center_angstrom"]),
            -float(record["n_h_center_angle_degrees"]),
            int(record["donor_id"]),
        )
    )

    shells: list[dict[str, Any]] = []
    for radius in radii:
        selected = [
            record
            for record in records
            if float(record["hydrogen_to_center_angstrom"]) <= radius
        ]
        inward = [record for record in selected if bool(record["inward_facing"])]
        shells.append(
            {
                "radius_angstrom": radius,
                "nh_count": len(selected),
                "inward_nh_count": len(inward),
                "component_counts": dict(
                    sorted(Counter(str(record.get("component")) for record in selected).items())
                ),
                "inward_component_counts": dict(
                    sorted(Counter(str(record.get("component")) for record in inward).items())
                ),
                "ranking_score": float(
                    sum(float(record["ranking_weight"]) for record in selected)
                ),
            }
        )
    return {
        "definition": (
            "covalent N-H groups around a geometric/electronic center; association proxy only, "
            "not proof of a classical hydrogen bond"
        ),
        "point_angstrom": point.tolist(),
        "nh_donor_count_total": len(records),
        "inward_angle_minimum_degrees": float(inward_angle_minimum_degrees),
        "distance_decay_angstrom": float(distance_decay_angstrom),
        "shells": shells,
        "nearest_donors": records[: min(12, len(records))],
    }


def select_nh_contrast_sites(
    sites: Sequence[Mapping[str, Any]],
    donors: Sequence[Mapping[str, Any]],
    positions_angstrom: NDArray[np.float64],
    cell_angstrom: NDArray[np.float64],
    *,
    shell_radii_angstrom: Sequence[float],
    selection_radius_angstrom: float,
    inward_angle_minimum_degrees: float,
    distance_decay_angstrom: float,
    minimum_separation_angstrom: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Choose separated N-H-facing and control sites from existing voids."""

    if len(sites) < 2:
        raise ValueError("at least two void sites are required")
    if minimum_separation_angstrom <= 0:
        raise ValueError("minimum site separation must be positive")
    radii = sorted({float(value) for value in shell_radii_angstrom})
    if selection_radius_angstrom not in radii:
        raise ValueError("N-H selection radius must be one of the shell radii")

    scored: list[dict[str, Any]] = []
    for site in sites:
        environment = nh_environment(
            site["cartesian_angstrom"],
            donors,
            positions_angstrom,
            cell_angstrom,
            shell_radii_angstrom=radii,
            inward_angle_minimum_degrees=inward_angle_minimum_degrees,
            distance_decay_angstrom=distance_decay_angstrom,
        )
        shell = next(
            item
            for item in environment["shells"]
            if float(item["radius_angstrom"]) == float(selection_radius_angstrom)
        )
        scored.append(
            {
                "site": dict(site),
                "nh_environment": environment,
                "selection_inward_nh_count": int(shell["inward_nh_count"]),
                "selection_nh_count": int(shell["nh_count"]),
                "selection_ranking_score": float(shell["ranking_score"]),
            }
        )
    scored.sort(key=lambda record: int(record["site"]["rank"]))
    facing = min(
        scored,
        key=lambda record: (
            -int(record["selection_inward_nh_count"]),
            -float(record["selection_ranking_score"]),
            -int(record["selection_nh_count"]),
            int(record["site"]["rank"]),
        ),
    )
    eligible = [
        record
        for record in scored
        if record is not facing
        and float(
            np.linalg.norm(
                minimum_image_vector(
                    record["site"]["cartesian_angstrom"],
                    facing["site"]["cartesian_angstrom"],
                    cell_angstrom,
                )
            )
        )
        >= minimum_separation_angstrom
    ]
    if not eligible:
        raise ValueError("no N-H control site satisfies the minimum separation")
    control = min(
        eligible,
        key=lambda record: (
            int(record["selection_inward_nh_count"]),
            float(record["selection_ranking_score"]),
            int(record["selection_nh_count"]),
            int(record["site"]["rank"]),
        ),
    )
    degenerate = (
        int(facing["selection_inward_nh_count"])
        == int(control["selection_inward_nh_count"])
        and abs(
            float(facing["selection_ranking_score"])
            - float(control["selection_ranking_score"])
        )
        <= 1.0e-12
    )

    def finalize(role: str, record: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "seed_role": role,
            "site": record["site"],
            "nh_environment": record["nh_environment"],
            "selection_radius_angstrom": float(selection_radius_angstrom),
            "selection_inward_nh_count": int(record["selection_inward_nh_count"]),
            "selection_nh_count": int(record["selection_nh_count"]),
            "selection_ranking_score": float(record["selection_ranking_score"]),
            "nh_degenerate_reference": degenerate,
        }

    return finalize("nh_facing", facing), finalize("nh_control", control)
