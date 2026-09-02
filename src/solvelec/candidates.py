from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from .trajectory import (
    VDW_RADII_ANGSTROM,
    fractional_grid,
    infer_element,
    minimum_image_vectors,
    surface_clearance,
)

_CANDIDATE_ID = re.compile(r"^[a-z][a-z0-9_-]*$")


def read_xyz(path: str | Path) -> tuple[list[str], np.ndarray, str]:
    """Read the small XYZ subset used for immutable Stage-A snapshots."""

    xyz_path = Path(path)
    lines = xyz_path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"invalid XYZ file: {xyz_path}")
    try:
        atom_count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"invalid XYZ atom count: {xyz_path}") from exc
    if atom_count <= 0 or len(lines) != atom_count + 2:
        raise ValueError(
            f"XYZ line count does not match atom count in {xyz_path}: "
            f"expected {atom_count + 2}, found {len(lines)}"
        )
    elements: list[str] = []
    positions: list[list[float]] = []
    for line_number, line in enumerate(lines[2:], start=3):
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"invalid XYZ atom record at {xyz_path}:{line_number}")
        label = fields[0].strip()
        element = "GH" if label.casefold() == "gh" else infer_element(label)
        try:
            position = [float(value) for value in fields[1:4]]
        except ValueError as exc:
            raise ValueError(f"invalid XYZ coordinates at {xyz_path}:{line_number}") from exc
        if not np.all(np.isfinite(position)):
            raise ValueError(f"non-finite XYZ coordinates at {xyz_path}:{line_number}")
        elements.append(element)
        positions.append(position)
    return elements, np.asarray(positions, dtype=float), lines[1]


def write_xyz(
    path: str | Path,
    elements: Sequence[str],
    positions: np.ndarray,
    comment: str,
) -> None:
    coordinates = np.asarray(positions, dtype=float)
    if coordinates.shape != (len(elements), 3) or not np.all(np.isfinite(coordinates)):
        raise ValueError("XYZ elements and coordinates are inconsistent")
    lines = [str(len(elements)), comment]
    for element, position in zip(elements, coordinates, strict=True):
        label = "Li" if element.upper() == "LI" else element
        lines.append(f"{label:<2s} {position[0]: .10f} {position[1]: .10f} {position[2]: .10f}")
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(output)


def _periodic_distance(first: np.ndarray, second: np.ndarray, cell: np.ndarray) -> float:
    vector = minimum_image_vectors(first[None, :], second[None, :], cell)[0, 0]
    return float(np.linalg.norm(vector))


def _refine_site(
    fractional: np.ndarray,
    atom_positions: np.ndarray,
    atom_radii: np.ndarray,
    cell: np.ndarray,
    *,
    points_per_axis: int,
    refinement_levels: int,
) -> tuple[np.ndarray, float]:
    offsets = np.asarray(
        [(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)],
        dtype=float,
    )
    best = np.asarray(fractional, dtype=float)
    clearance = float(surface_clearance(best[None, :], atom_positions, atom_radii, cell)[0])
    step = 0.5 / points_per_axis
    for _ in range(refinement_levels):
        trial = (best[None, :] + offsets * step) % 1.0
        values = surface_clearance(trial, atom_positions, atom_radii, cell)
        index = int(np.argmax(values))
        best = trial[index]
        clearance = float(values[index])
        step /= 3.0
    return best, clearance


def ranked_void_sites(
    atom_positions: np.ndarray,
    atom_radii: np.ndarray,
    cell: np.ndarray,
    *,
    count: int,
    points_per_axis: int,
    refinement_levels: int,
    minimum_separation_angstrom: float,
    minimum_clearance_angstrom: float,
) -> list[dict[str, Any]]:
    """Return deterministic, non-overlapping local maxima of a void proxy."""

    if count < 2 or points_per_axis < 2 or refinement_levels < 0:
        raise ValueError("void search requires count/axis >= 2 and non-negative refinement")
    if minimum_separation_angstrom <= 0 or minimum_clearance_angstrom < 0:
        raise ValueError("void separation must be positive and clearance non-negative")
    matrix = np.asarray(cell, dtype=float)
    grid = fractional_grid(points_per_axis)
    clearances = surface_clearance(grid, atom_positions, atom_radii, matrix)
    order = np.lexsort((grid[:, 2], grid[:, 1], grid[:, 0], -clearances))
    seed_count = min(len(grid), max(count * 8, count))
    refined: list[tuple[np.ndarray, float]] = []
    for index in order:
        if clearances[index] < minimum_clearance_angstrom:
            break
        cartesian = grid[index] @ matrix
        if any(
            _periodic_distance(cartesian, previous_fractional @ matrix, matrix)
            < minimum_separation_angstrom * 0.75
            for previous_fractional, _ in refined
        ):
            continue
        refined.append(
            _refine_site(
                grid[index],
                atom_positions,
                atom_radii,
                matrix,
                points_per_axis=points_per_axis,
                refinement_levels=refinement_levels,
            )
        )
        if len(refined) >= seed_count:
            break
    refined.sort(key=lambda item: (-item[1], *item[0].tolist()))
    accepted: list[tuple[np.ndarray, float]] = []
    for fractional, clearance in refined:
        cartesian = fractional @ matrix
        if clearance < minimum_clearance_angstrom:
            continue
        if any(
            _periodic_distance(cartesian, previous @ matrix, matrix)
            < minimum_separation_angstrom
            for previous, _ in accepted
        ):
            continue
        accepted.append((fractional, clearance))
        if len(accepted) == count:
            break
    if len(accepted) < count:
        raise ValueError(
            f"void search found {len(accepted)} sites but {count} are required; "
            "lower the configured site count/separation only after inspection"
        )
    return [
        {
            "rank": rank,
            "radius_angstrom": float(clearance),
            "fractional": fractional.tolist(),
            "cartesian_angstrom": (fractional @ matrix).tolist(),
        }
        for rank, (fractional, clearance) in enumerate(accepted, start=1)
    ]


def select_candidate_pairs(
    sites: Sequence[Mapping[str, Any]],
    cell: np.ndarray,
    definitions: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Map configured Li--cavity target distances to distinct void-site pairs."""

    if len(sites) < 2:
        raise ValueError("at least two void sites are required")
    matrix = np.asarray(cell, dtype=float)
    used_pairs: set[tuple[int, int]] = set()
    selected: list[dict[str, Any]] = []
    for definition in definitions:
        candidate_id = str(definition["id"])
        if not _CANDIDATE_ID.fullmatch(candidate_id):
            raise ValueError(f"invalid Stage-B candidate id {candidate_id!r}")
        target = float(definition["target_li_cavity_distance_angstrom"])
        tolerance = float(definition["tolerance_angstrom"])
        if target <= 0 or tolerance <= 0:
            raise ValueError(f"candidate {candidate_id} target/tolerance must be positive")
        options: list[tuple[tuple[float, float, float, int, int], int, int, float]] = []
        for first in range(len(sites)):
            for second in range(first + 1, len(sites)):
                if (first, second) in used_pairs:
                    continue
                first_position = np.asarray(sites[first]["cartesian_angstrom"], dtype=float)
                second_position = np.asarray(sites[second]["cartesian_angstrom"], dtype=float)
                distance = _periodic_distance(first_position, second_position, matrix)
                error = abs(distance - target)
                if error > tolerance:
                    continue
                first_clearance = float(sites[first]["radius_angstrom"])
                second_clearance = float(sites[second]["radius_angstrom"])
                score = (
                    error,
                    -min(first_clearance, second_clearance),
                    -(first_clearance + second_clearance),
                    first,
                    second,
                )
                options.append((score, first, second, distance))
        if not options:
            raise ValueError(
                f"no distinct void pair satisfies candidate {candidate_id!r}: "
                f"target {target:.3f} +/- {tolerance:.3f} angstrom"
            )
        _, first, second, distance = min(options)
        used_pairs.add((first, second))
        # The roomier site receives Li. This is only a deterministic insertion
        # seed; neither site is interpreted as an observed electron position.
        if float(sites[second]["radius_angstrom"]) > float(
            sites[first]["radius_angstrom"]
        ):
            first, second = second, first
        selected.append(
            {
                "id": candidate_id,
                "target_li_cavity_distance_angstrom": target,
                "tolerance_angstrom": tolerance,
                "achieved_li_cavity_distance_angstrom": distance,
                "li_site": dict(sites[first]),
                "cavity_basis_site": dict(sites[second]),
            }
        )
    return selected


def heavy_atom_geometry(
    elements: Sequence[str], positions: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    heavy_positions: list[np.ndarray] = []
    radii: list[float] = []
    for element, position in zip(elements, np.asarray(positions, dtype=float), strict=True):
        normalized = infer_element(element)
        if normalized == "H":
            continue
        heavy_positions.append(position)
        radii.append(VDW_RADII_ANGSTROM[normalized])
    if not heavy_positions:
        raise ValueError("candidate construction requires at least one heavy atom")
    return np.asarray(heavy_positions, dtype=float), np.asarray(radii, dtype=float)
