from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .pbc import periodic_radius_of_gyration, periodic_weighted_centroid


@dataclass(frozen=True)
class CubeData:
    comments: tuple[str, str]
    origin: NDArray[np.float64]
    axes: NDArray[np.float64]
    atoms: NDArray[np.float64]
    values: NDArray[np.float64]
    coordinate_unit: str

    @property
    def shape(self) -> tuple[int, int, int]:
        shape = self.values.shape
        if len(shape) != 3:
            raise ValueError(f"Cube density must be three-dimensional, got shape {shape}")
        return (int(shape[0]), int(shape[1]), int(shape[2]))

    @property
    def box(self) -> NDArray[np.float64]:
        counts = np.asarray(self.shape, dtype=float)
        return self.axes * counts[:, None]

    @property
    def voxel_volume(self) -> float:
        return float(abs(np.linalg.det(self.axes)))

    def grid_positions(self) -> NDArray[np.float64]:
        indices = np.indices(self.shape, dtype=float).reshape(3, -1).T
        return self.origin + indices @ self.axes


@dataclass(frozen=True)
class SpinDensityMetrics:
    electron_count: float
    signed_integral: float
    centroid: tuple[float, float, float]
    radius: float
    inverse_participation_ratio: float
    positive_voxels: int

    def as_dict(self) -> dict[str, float | int | list[float]]:
        return {
            "electron_count": self.electron_count,
            "signed_integral": self.signed_integral,
            "centroid": list(self.centroid),
            "radius": self.radius,
            "inverse_participation_ratio": self.inverse_participation_ratio,
            "positive_voxels": self.positive_voxels,
        }


def read_cube(path: str | Path) -> CubeData:
    cube_path = Path(path)
    with cube_path.open("r", encoding="utf-8") as handle:
        comments = (handle.readline().rstrip(), handle.readline().rstrip())
        header = handle.readline().split()
        if len(header) < 4:
            raise ValueError(f"Invalid cube atom/origin line in {cube_path}")
        atom_count = abs(int(header[0]))
        origin = np.asarray([float(value) for value in header[1:4]], dtype=float)
        counts: list[int] = []
        axes: list[list[float]] = []
        signed_counts: list[int] = []
        for _ in range(3):
            fields = handle.readline().split()
            if len(fields) < 4:
                raise ValueError(f"Invalid cube grid line in {cube_path}")
            signed = int(fields[0])
            signed_counts.append(signed)
            counts.append(abs(signed))
            axes.append([float(value) for value in fields[1:4]])
        atoms = []
        for _ in range(atom_count):
            fields = handle.readline().split()
            if len(fields) < 5:
                raise ValueError(f"Invalid cube atom record in {cube_path}")
            atoms.append([float(value) for value in fields[:5]])
        values = np.fromstring(handle.read(), sep=" ", dtype=float)
    expected = int(np.prod(counts))
    if values.size != expected:
        raise ValueError(f"Cube contains {values.size} values; expected {expected}")
    unit = "angstrom" if all(value < 0 for value in signed_counts) else "bohr"
    return CubeData(
        comments=comments,
        origin=origin,
        axes=np.asarray(axes, dtype=float),
        atoms=np.asarray(atoms, dtype=float).reshape((-1, 5)),
        values=values.reshape(tuple(counts)),
        coordinate_unit=unit,
    )


def analyze_spin_density(cube: CubeData, clip_negative: bool = True) -> SpinDensityMetrics:
    density = cube.values.reshape(-1)
    signed_integral = float(density.sum() * cube.voxel_volume)
    weights_density = np.clip(density, 0.0, None) if clip_negative else np.abs(density)
    weights = weights_density * cube.voxel_volume
    electron_count = float(weights.sum())
    if electron_count <= 0:
        raise ValueError("Spin-density cube has no positive integrated density")
    positions = cube.grid_positions()
    centroid = periodic_weighted_centroid(positions, weights, cube.box, cube.origin)
    radius = periodic_radius_of_gyration(positions, weights, centroid, cube.box)
    ipr = float(np.sum(weights_density**2) * cube.voxel_volume / electron_count**2)
    return SpinDensityMetrics(
        electron_count=electron_count,
        signed_integral=signed_integral,
        centroid=(float(centroid[0]), float(centroid[1]), float(centroid[2])),
        radius=radius,
        inverse_participation_ratio=ipr,
        positive_voxels=int(np.count_nonzero(weights_density > 0)),
    )
