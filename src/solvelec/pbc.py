from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def as_box_matrix(box: ArrayLike) -> NDArray[np.float64]:
    matrix = np.asarray(box, dtype=float)
    if matrix.shape == (3,):
        matrix = np.diag(matrix)
    if matrix.shape != (3, 3):
        raise ValueError("box must be three lengths or a 3x3 matrix of cell vectors")
    if abs(np.linalg.det(matrix)) < 1.0e-12:
        raise ValueError("box matrix is singular")
    return matrix


def minimum_image_displacement(
    point: ArrayLike, reference: ArrayLike, box: ArrayLike
) -> NDArray[np.float64]:
    matrix = as_box_matrix(box)
    delta = np.asarray(point, dtype=float) - np.asarray(reference, dtype=float)
    fractional = delta @ np.linalg.inv(matrix)
    fractional -= np.floor(fractional + 0.5)
    return fractional @ matrix


def periodic_weighted_centroid(
    positions: ArrayLike, weights: ArrayLike, box: ArrayLike, origin: ArrayLike | None = None
) -> NDArray[np.float64]:
    coords = np.asarray(positions, dtype=float)
    weight = np.asarray(weights, dtype=float)
    matrix = as_box_matrix(box)
    cell_origin = np.zeros(3) if origin is None else np.asarray(origin, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("positions must have shape (n, 3)")
    if weight.shape != (coords.shape[0],):
        raise ValueError("weights must have one value per position")
    total = float(weight.sum())
    if total <= 0:
        raise ValueError("weights must sum to a positive value")
    fractional = (coords - cell_origin) @ np.linalg.inv(matrix)
    fractional %= 1.0
    centroid_fractional = np.empty(3)
    for axis in range(3):
        angles = 2.0 * np.pi * fractional[:, axis]
        sine = float(np.dot(weight, np.sin(angles)))
        cosine = float(np.dot(weight, np.cos(angles)))
        if abs(sine) + abs(cosine) < 1.0e-14:
            raise ValueError("periodic centroid is undefined for a symmetric distribution")
        angle = np.arctan2(sine, cosine)
        centroid_fractional[axis] = (angle / (2.0 * np.pi)) % 1.0
    return cell_origin + centroid_fractional @ matrix


def periodic_radius_of_gyration(
    positions: ArrayLike,
    weights: ArrayLike,
    centroid: ArrayLike,
    box: ArrayLike,
) -> float:
    coords = np.asarray(positions, dtype=float)
    weight = np.asarray(weights, dtype=float)
    total = float(weight.sum())
    if total <= 0:
        raise ValueError("weights must sum to a positive value")
    displacements = np.vstack(
        [minimum_image_displacement(point, centroid, box) for point in coords]
    )
    return float(
        np.sqrt(np.dot(weight, np.einsum("ij,ij->i", displacements, displacements)) / total)
    )
