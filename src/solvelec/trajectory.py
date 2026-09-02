from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np

VDW_RADII_ANGSTROM = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "LI": 1.82,
}


def infer_element(atom_name: str, explicit_element: str | None = None) -> str:
    """Return a conservative element symbol from topology metadata."""

    if explicit_element:
        normalized = explicit_element.strip().upper()
        if normalized in VDW_RADII_ANGSTROM:
            return normalized
    letters = "".join(character for character in atom_name if character.isalpha()).upper()
    if letters.startswith("LI"):
        return "LI"
    if letters and letters[0] in VDW_RADII_ANGSTROM:
        return letters[0]
    raise ValueError(f"cannot infer supported element from atom name {atom_name!r}")


def cell_matrix(dimensions: Sequence[float]) -> np.ndarray:
    """Convert ``a,b,c,alpha,beta,gamma`` into row cell vectors in Angstrom."""

    values = np.asarray(dimensions, dtype=float)
    if values.shape != (6,) or not np.all(np.isfinite(values)):
        raise ValueError("cell dimensions must contain six finite values")
    a, b, c, alpha_deg, beta_deg, gamma_deg = values
    if min(a, b, c) <= 0 or not all(0 < angle < 180 for angle in values[3:]):
        raise ValueError("cell lengths and angles must be physically valid")
    alpha, beta, gamma = np.deg2rad([alpha_deg, beta_deg, gamma_deg])
    sin_gamma = float(np.sin(gamma))
    if abs(sin_gamma) < 1.0e-8:
        raise ValueError("cell gamma angle is singular")
    vector_a = np.array([a, 0.0, 0.0])
    vector_b = np.array([b * np.cos(gamma), b * sin_gamma, 0.0])
    vector_c_x = c * np.cos(beta)
    vector_c_y = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / sin_gamma
    vector_c_z_sq = c**2 - vector_c_x**2 - vector_c_y**2
    if vector_c_z_sq <= 0:
        raise ValueError("cell vectors have a non-positive volume")
    matrix = np.vstack(
        [vector_a, vector_b, np.array([vector_c_x, vector_c_y, np.sqrt(vector_c_z_sq)])]
    )
    if abs(float(np.linalg.det(matrix))) < 1.0e-8:
        raise ValueError("cell matrix is singular")
    return matrix


def fractional_grid(points_per_axis: int) -> np.ndarray:
    if points_per_axis < 2:
        raise ValueError("at least two grid points per axis are required")
    axis = (np.arange(points_per_axis, dtype=float) + 0.5) / points_per_axis
    mesh = np.meshgrid(axis, axis, axis, indexing="ij")
    return np.column_stack([component.ravel() for component in mesh])


def minimum_image_vectors(origins: np.ndarray, targets: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Return pairwise target-origin vectors under periodic boundary conditions."""

    origin_array = np.asarray(origins, dtype=float)
    target_array = np.asarray(targets, dtype=float)
    matrix = np.asarray(cell, dtype=float)
    if origin_array.ndim != 2 or origin_array.shape[1] != 3:
        raise ValueError("origins must be shaped (n, 3)")
    if target_array.ndim != 2 or target_array.shape[1] != 3:
        raise ValueError("targets must be shaped (m, 3)")
    inverse = np.linalg.inv(matrix)
    fractional = (target_array[None, :, :] - origin_array[:, None, :]) @ inverse
    fractional -= np.rint(fractional)
    return fractional @ matrix


def pair_distances(
    first: np.ndarray,
    second: np.ndarray,
    cell: np.ndarray,
    *,
    same_group: bool = False,
    first_residue_ids: np.ndarray | None = None,
    second_residue_ids: np.ndarray | None = None,
    exclude_same_residue: bool = False,
) -> np.ndarray:
    vectors = minimum_image_vectors(first, second, cell)
    distances = np.linalg.norm(vectors, axis=2)
    mask = np.ones(distances.shape, dtype=bool)
    if same_group:
        if distances.shape[0] != distances.shape[1]:
            raise ValueError("same-group distances require equal-sized groups")
        mask &= np.triu(np.ones(distances.shape, dtype=bool), k=1)
    if exclude_same_residue:
        if first_residue_ids is None or second_residue_ids is None:
            raise ValueError("residue ids are required for same-residue exclusion")
        mask &= np.asarray(first_residue_ids)[:, None] != np.asarray(second_residue_ids)[None, :]
    return distances[mask]


def surface_clearance(
    fractional_points: np.ndarray,
    atom_positions: np.ndarray,
    atom_radii: np.ndarray,
    cell: np.ndarray,
    *,
    chunk_size: int = 256,
) -> np.ndarray:
    """Nearest heavy-atom van-der-Waals surface distance for fractional points."""

    points = np.asarray(fractional_points, dtype=float)
    atoms = np.asarray(atom_positions, dtype=float)
    radii = np.asarray(atom_radii, dtype=float)
    matrix = np.asarray(cell, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) == 0:
        raise ValueError("fractional points must be a non-empty (n, 3) array")
    if atoms.ndim != 2 or atoms.shape[1] != 3 or len(atoms) == 0:
        raise ValueError("atom positions must be a non-empty (m, 3) array")
    if radii.shape != (len(atoms),) or np.any(radii <= 0):
        raise ValueError("one positive radius is required per atom")
    inverse = np.linalg.inv(matrix)
    atom_fractional = atoms @ inverse
    values = np.empty(len(points), dtype=float)
    for start in range(0, len(points), chunk_size):
        stop = min(start + chunk_size, len(points))
        difference = atom_fractional[None, :, :] - points[start:stop, None, :]
        difference -= np.rint(difference)
        distances = np.linalg.norm(difference @ matrix, axis=2)
        values[start:stop] = np.min(distances - radii[None, :], axis=1)
    return values


def largest_void_proxy(
    atom_positions: np.ndarray,
    atom_radii: np.ndarray,
    cell: np.ndarray,
    *,
    points_per_axis: int,
    refinement_levels: int = 2,
) -> dict[str, Any]:
    """Find a deterministic grid/refinement proxy for the largest interstitial void."""

    if refinement_levels < 0:
        raise ValueError("refinement_levels must be non-negative")
    matrix = np.asarray(cell, dtype=float)
    candidates = fractional_grid(points_per_axis)
    clearances = surface_clearance(candidates, atom_positions, atom_radii, matrix)
    best_index = int(np.argmax(clearances))
    best_fractional = candidates[best_index]
    best_clearance = float(clearances[best_index])
    step = 0.5 / points_per_axis
    offsets_base = np.asarray(list(_offset_triplets()), dtype=float)
    for _ in range(refinement_levels):
        candidates = (best_fractional[None, :] + offsets_base * step) % 1.0
        clearances = surface_clearance(candidates, atom_positions, atom_radii, matrix)
        best_index = int(np.argmax(clearances))
        best_fractional = candidates[best_index]
        best_clearance = float(clearances[best_index])
        step /= 3.0
    best_cartesian = best_fractional @ matrix
    return {
        "radius_angstrom": best_clearance,
        "fractional": best_fractional.tolist(),
        "cartesian_angstrom": best_cartesian.tolist(),
    }


def _offset_triplets() -> Iterable[tuple[int, int, int]]:
    for first in (-1, 0, 1):
        for second in (-1, 0, 1):
            for third in (-1, 0, 1):
                yield first, second, third


def autocorrelation_summary(times_ps: Sequence[float], values: Sequence[float]) -> dict[str, float]:
    """Estimate statistical inefficiency using the initial-positive ACF sequence."""

    times = np.asarray(times_ps, dtype=float)
    series = np.asarray(values, dtype=float)
    if times.ndim != 1 or series.ndim != 1 or len(times) != len(series) or len(times) < 4:
        raise ValueError("times and values must be same-length vectors with at least four points")
    if not np.all(np.isfinite(times)) or not np.all(np.diff(times) > 0):
        raise ValueError("times must be finite and strictly increasing")
    if not np.all(np.isfinite(series)):
        raise ValueError("values must be finite")
    centered = series - np.mean(series)
    variance = float(np.dot(centered, centered) / len(centered))
    spacing = float(np.median(np.diff(times)))
    if variance <= np.finfo(float).eps:
        return {
            "statistical_inefficiency": 1.0,
            "integrated_autocorrelation_time_ps": 0.5 * spacing,
            "effective_sample_size": float(len(series)),
        }
    padded = 1 << (2 * len(series) - 1).bit_length()
    transform = np.fft.rfft(centered, n=padded)
    autocovariance = np.fft.irfft(transform * np.conjugate(transform), n=padded)[: len(series)]
    autocovariance /= np.arange(len(series), 0, -1)
    correlation = autocovariance / autocovariance[0]
    positive = correlation[1:]
    stop = np.flatnonzero(positive <= 0)
    if len(stop):
        positive = positive[: int(stop[0])]
    inefficiency = max(1.0, 1.0 + 2.0 * float(np.sum(positive)))
    return {
        "statistical_inefficiency": inefficiency,
        "integrated_autocorrelation_time_ps": 0.5 * inefficiency * spacing,
        "effective_sample_size": float(len(series) / inefficiency),
    }


def select_representative_indices(
    records: Sequence[dict[str, float]],
    feature_names: Sequence[str],
    *,
    count: int,
    minimum_time_ps: float,
    minimum_separation_ps: float,
) -> list[int]:
    """Select deterministic medoid/farthest-point representatives from an equilibrated pool."""

    if count <= 0 or minimum_separation_ps < 0:
        raise ValueError("count must be positive and separation non-negative")
    eligible = [index for index, row in enumerate(records) if row["time_ps"] >= minimum_time_ps]
    if len(eligible) < count:
        raise ValueError("not enough equilibrated records for the requested snapshot count")
    features = np.asarray(
        [[float(records[index][name]) for name in feature_names] for index in eligible],
        dtype=float,
    )
    if not np.all(np.isfinite(features)):
        raise ValueError("selection features must be finite")
    center = np.median(features, axis=0)
    scale = np.median(np.abs(features - center), axis=0) * 1.4826
    standard_deviation = np.std(features, axis=0)
    scale = np.where(scale > 1.0e-12, scale, standard_deviation)
    scale = np.where(scale > 1.0e-12, scale, 1.0)
    standardized = (features - center) / scale
    center_distance = np.linalg.norm(standardized, axis=1)
    first_local = min(
        range(len(eligible)),
        key=lambda local: (
            float(center_distance[local]),
            float(records[eligible[local]]["time_ps"]),
            eligible[local],
        ),
    )
    chosen_local = [first_local]
    while len(chosen_local) < count:
        candidates: list[tuple[float, float, int]] = []
        for local, global_index in enumerate(eligible):
            if local in chosen_local:
                continue
            time = float(records[global_index]["time_ps"])
            if any(
                abs(time - float(records[eligible[selected]]["time_ps"])) < minimum_separation_ps
                for selected in chosen_local
            ):
                continue
            diversity = min(
                float(np.linalg.norm(standardized[local] - standardized[selected]))
                for selected in chosen_local
            )
            candidates.append((diversity, -time, local))
        if not candidates:
            raise ValueError("decorrelation separation leaves too few snapshot candidates")
        chosen_local.append(max(candidates)[2])
    return [eligible[local] for local in sorted(chosen_local, key=lambda item: eligible[item])]
