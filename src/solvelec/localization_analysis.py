from __future__ import annotations

from collections import Counter, deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .cube import CubeData, analyze_spin_density
from .trajectory import VDW_RADII_ANGSTROM, minimum_image_vectors

BOHR_TO_ANGSTROM = 0.529177210903
COVALENT_RADII_ANGSTROM = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
}


@dataclass(frozen=True)
class MolecularTopology:
    molecule_ids: NDArray[np.int64]
    molecules: tuple[dict[str, Any], ...]
    problems: tuple[str, ...]


@dataclass(frozen=True)
class GeometricPartition:
    owner: NDArray[np.int32]
    cavity_indices: NDArray[np.int64]
    atom_radii_angstrom: NDArray[np.float64]


def read_cp2k_cell(path: str | Path) -> NDArray[np.float64]:
    vectors: dict[str, list[float]] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) != 4 or fields[0].upper() not in {"A", "B", "C"}:
            continue
        vectors[fields[0].upper()] = [float(value) for value in fields[1:4]]
    if set(vectors) != {"A", "B", "C"}:
        raise ValueError(f"cell include {path} must define A, B, and C vectors")
    matrix = np.asarray([vectors[axis] for axis in ("A", "B", "C")], dtype=float)
    if not np.all(np.isfinite(matrix)) or float(np.linalg.det(matrix)) <= 0:
        raise ValueError(f"cell include {path} has an invalid cell matrix")
    return matrix


def length_scale_to_angstrom(cube: CubeData) -> float:
    if cube.coordinate_unit == "bohr":
        return BOHR_TO_ANGSTROM
    if cube.coordinate_unit == "angstrom":
        return 1.0
    raise ValueError(f"unsupported cube coordinate unit {cube.coordinate_unit!r}")


def cube_origin_angstrom(cube: CubeData) -> NDArray[np.float64]:
    return np.asarray(cube.origin, dtype=float) * length_scale_to_angstrom(cube)


def cube_axes_angstrom(cube: CubeData) -> NDArray[np.float64]:
    return np.asarray(cube.axes, dtype=float) * length_scale_to_angstrom(cube)


def cube_cell_angstrom(cube: CubeData) -> NDArray[np.float64]:
    counts = np.asarray(cube.shape, dtype=float)
    return cube_axes_angstrom(cube) * counts[:, None]


def cube_atom_positions_angstrom(cube: CubeData) -> NDArray[np.float64]:
    return np.asarray(cube.atoms[:, 2:5], dtype=float) * length_scale_to_angstrom(cube)


def cube_compatibility_problems(
    first: CubeData,
    second: CubeData,
    *,
    tolerance_angstrom: float,
) -> list[str]:
    problems: list[str] = []
    if first.shape != second.shape:
        problems.append(f"cube shapes differ: {first.shape} versus {second.shape}")
        return problems
    if not np.allclose(
        cube_origin_angstrom(first),
        cube_origin_angstrom(second),
        atol=tolerance_angstrom,
        rtol=0.0,
    ):
        problems.append("cube origins differ")
    if not np.allclose(
        cube_axes_angstrom(first),
        cube_axes_angstrom(second),
        atol=tolerance_angstrom,
        rtol=0.0,
    ):
        problems.append("cube grid vectors differ")
    if len(first.atoms) != len(second.atoms):
        problems.append(f"cube atom counts differ: {len(first.atoms)} versus {len(second.atoms)}")
    elif len(first.atoms):
        if not np.array_equal(first.atoms[:, 0].astype(int), second.atoms[:, 0].astype(int)):
            problems.append("cube atomic-number sequences differ")
        if not np.allclose(
            cube_atom_positions_angstrom(first),
            cube_atom_positions_angstrom(second),
            atol=tolerance_angstrom,
            rtol=0.0,
        ):
            problems.append("cube atom positions differ")
    return problems


def candidate_cube_problems(
    cube: CubeData,
    elements: Sequence[str],
    positions_angstrom: NDArray[np.float64],
    cell_angstrom: NDArray[np.float64],
    *,
    tolerance_angstrom: float,
) -> list[str]:
    problems: list[str] = []
    cube_cell = cube_cell_angstrom(cube)
    if not np.allclose(cube_cell, cell_angstrom, atol=tolerance_angstrom, rtol=0.0):
        maximum = float(np.max(np.abs(cube_cell - cell_angstrom)))
        problems.append(f"cube and candidate cells differ by up to {maximum:.6g} angstrom")
    real_indices = [index for index, element in enumerate(elements) if element.upper() != "GH"]
    candidates = np.asarray(positions_angstrom, dtype=float)
    if len(cube.atoms) == len(elements):
        reference = candidates
    elif len(cube.atoms) == len(real_indices):
        reference = candidates[real_indices]
    else:
        problems.append(
            "cube atom count does not match candidate with or without its ghost center: "
            f"{len(cube.atoms)} versus {len(elements)}/{len(real_indices)}"
        )
        return problems
    cube_positions = cube_atom_positions_angstrom(cube)
    fractional = (reference - cube_positions) @ np.linalg.inv(cell_angstrom)
    fractional -= np.rint(fractional)
    paired = fractional @ cell_angstrom
    maximum = float(np.max(np.linalg.norm(paired, axis=1))) if len(reference) else 0.0
    if maximum > tolerance_angstrom:
        problems.append(f"cube and candidate atom positions differ by up to {maximum:.6g} angstrom")
    return problems


def _formula(elements: Sequence[str]) -> dict[str, int]:
    return dict(sorted(Counter(element.upper() for element in elements).items()))


def infer_molecular_topology(
    elements: Sequence[str],
    positions_angstrom: NDArray[np.float64],
    cell_angstrom: NDArray[np.float64],
    component_counts: Mapping[str, int],
    component_definitions: Mapping[str, Mapping[str, Any]],
    *,
    bond_scale: float,
) -> MolecularTopology:
    if bond_scale <= 1.0:
        raise ValueError("bond scale must be greater than one")
    normalized = [element.upper() for element in elements]
    solvent_indices = [
        index for index, element in enumerate(normalized) if element not in {"LI", "GH"}
    ]
    if not solvent_indices:
        raise ValueError("candidate contains no solvent atoms")
    unsupported = sorted(
        {normalized[index] for index in solvent_indices} - set(COVALENT_RADII_ANGSTROM)
    )
    if unsupported:
        raise ValueError(f"unsupported elements for bond inference: {unsupported}")
    solvent_positions = np.asarray(positions_angstrom, dtype=float)[solvent_indices]
    vectors = minimum_image_vectors(solvent_positions, solvent_positions, cell_angstrom)
    distances = np.linalg.norm(vectors, axis=2)
    radii = np.asarray(
        [COVALENT_RADII_ANGSTROM[normalized[index]] for index in solvent_indices],
        dtype=float,
    )
    cutoffs = bond_scale * (radii[:, None] + radii[None, :])
    bonded = (distances > 0.35) & (distances <= cutoffs)
    np.fill_diagonal(bonded, False)

    unseen = set(range(len(solvent_indices)))
    components: list[list[int]] = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        queue: deque[int] = deque([seed])
        local_component: list[int] = []
        while queue:
            current = queue.popleft()
            local_component.append(current)
            neighbors = np.flatnonzero(bonded[current])
            for neighbor in neighbors:
                value = int(neighbor)
                if value in unseen:
                    unseen.remove(value)
                    queue.append(value)
        components.append(sorted(solvent_indices[index] for index in local_component))
    components.sort(key=lambda values: values[0])

    formulas_to_names: dict[tuple[tuple[str, int], ...], list[str]] = {}
    for name in component_counts:
        definition = component_definitions.get(name, {})
        formula = definition.get("element_counts")
        if not isinstance(formula, dict) or not formula:
            raise ValueError(f"component {name!r} lacks element_counts")
        key = tuple(
            sorted((str(element).upper(), int(count)) for element, count in formula.items())
        )
        formulas_to_names.setdefault(key, []).append(str(name))

    molecule_ids = np.full(len(elements), -1, dtype=np.int64)
    molecules: list[dict[str, Any]] = []
    detected_counts: Counter[str] = Counter()
    problems: list[str] = []
    for molecule_id, atom_indices in enumerate(components):
        formula = _formula([normalized[index] for index in atom_indices])
        key = tuple(formula.items())
        names = formulas_to_names.get(key, [])
        component = names[0] if len(names) == 1 else None
        if component is None:
            problems.append(
                f"molecule {molecule_id + 1} formula {formula} does not uniquely match "
                f"configured components {sorted(component_counts)}"
            )
        else:
            detected_counts[component] += 1
        molecule_ids[atom_indices] = molecule_id
        molecules.append(
            {
                "molecule_id": molecule_id + 1,
                "component": component,
                "formula": formula,
                "atom_indices_1based": [index + 1 for index in atom_indices],
            }
        )
    expected_counts = {str(name): int(count) for name, count in component_counts.items()}
    if dict(detected_counts) != expected_counts:
        problems.append(
            f"detected molecule counts {dict(sorted(detected_counts.items()))} differ from "
            f"specification {dict(sorted(expected_counts.items()))}"
        )
    return MolecularTopology(molecule_ids, tuple(molecules), tuple(problems))


def _sphere_grid_indices(
    cube: CubeData,
    center_angstrom: NDArray[np.float64],
    radius_angstrom: float,
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    if radius_angstrom <= 0:
        raise ValueError("sphere radius must be positive")
    origin = cube_origin_angstrom(cube)
    axes = cube_axes_angstrom(cube)
    inverse_axes = np.linalg.inv(axes)
    center_index = (np.asarray(center_angstrom, dtype=float) - origin) @ inverse_axes
    half_widths = np.ceil(
        radius_angstrom * np.linalg.norm(inverse_axes, axis=0)
    ).astype(int) + 1
    center_integer = np.rint(center_index).astype(int)
    ranges = [
        np.arange(-half_width, half_width + 1, dtype=int) for half_width in half_widths
    ]
    mesh = np.meshgrid(*ranges, indexing="ij")
    offsets = np.column_stack([item.reshape(-1) for item in mesh])
    counts = np.asarray(cube.shape, dtype=int)
    indices = (center_integer[None, :] + offsets) % counts[None, :]
    flat = np.ravel_multi_index(indices.T, cube.shape)
    flat = np.unique(flat).astype(np.int64)
    indices = np.column_stack(np.unravel_index(flat, cube.shape))
    positions = origin + indices @ axes
    vectors = minimum_image_vectors(
        np.asarray(center_angstrom, dtype=float).reshape(1, 3), positions, cube_cell_angstrom(cube)
    )[0]
    distances = np.linalg.norm(vectors, axis=1)
    inside = distances <= radius_angstrom
    return flat[inside], distances[inside]


def build_geometric_partition(
    cube: CubeData,
    elements: Sequence[str],
    positions_angstrom: NDArray[np.float64],
    cavity_center_angstrom: NDArray[np.float64],
    *,
    vdw_scale: float,
    cavity_probe_radius_angstrom: float,
) -> GeometricPartition:
    if vdw_scale <= 0:
        raise ValueError("van-der-Waals region scale must be positive")
    owner = np.full(int(np.prod(cube.shape)), -1, dtype=np.int32)
    best_scaled_distance = np.full(owner.shape, np.inf, dtype=np.float32)
    atom_radii = np.zeros(len(elements), dtype=float)
    for atom_index, (element, center) in enumerate(
        zip(elements, np.asarray(positions_angstrom, dtype=float), strict=True)
    ):
        normalized = element.upper()
        if normalized == "GH":
            continue
        if normalized not in VDW_RADII_ANGSTROM:
            raise ValueError(f"unsupported element for geometric partition: {element!r}")
        radius = float(VDW_RADII_ANGSTROM[normalized]) * vdw_scale
        atom_radii[atom_index] = radius
        indices, distances = _sphere_grid_indices(cube, center, radius)
        scaled = distances / radius
        update = scaled < best_scaled_distance[indices]
        selected = indices[update]
        owner[selected] = atom_index
        best_scaled_distance[selected] = scaled[update].astype(np.float32)
    cavity_indices, _ = _sphere_grid_indices(
        cube, np.asarray(cavity_center_angstrom, dtype=float), cavity_probe_radius_angstrom
    )
    return GeometricPartition(owner, cavity_indices, atom_radii)


def integrate_partitioned_field(
    cube: CubeData,
    field: NDArray[np.float64],
    partition: GeometricPartition,
    molecule_ids: NDArray[np.int64],
) -> dict[str, Any]:
    values = np.asarray(field, dtype=float).reshape(-1)
    if values.size != partition.owner.size:
        raise ValueError("field and geometric partition sizes differ")
    if molecule_ids.ndim != 1:
        raise ValueError("molecule ids must be one-dimensional")
    voxel = cube.voxel_volume
    positive = np.clip(values, 0.0, None) * voxel
    negative = np.clip(-values, 0.0, None) * voxel
    signed = values * voxel
    owned = partition.owner >= 0
    owner_values = partition.owner[owned].astype(int)
    atom_count = len(molecule_ids)
    atom_positive = np.bincount(owner_values, weights=positive[owned], minlength=atom_count)
    atom_negative = np.bincount(owner_values, weights=negative[owned], minlength=atom_count)
    atom_signed = np.bincount(owner_values, weights=signed[owned], minlength=atom_count)
    molecule_count = int(np.max(molecule_ids)) + 1 if np.any(molecule_ids >= 0) else 0
    molecule_positive = np.zeros(molecule_count, dtype=float)
    molecule_negative = np.zeros(molecule_count, dtype=float)
    molecule_signed = np.zeros(molecule_count, dtype=float)
    for atom_index, molecule_id in enumerate(molecule_ids):
        if molecule_id < 0:
            continue
        molecule_positive[molecule_id] += atom_positive[atom_index]
        molecule_negative[molecule_id] += atom_negative[atom_index]
        molecule_signed[molecule_id] += atom_signed[atom_index]
    positive_total = float(positive.sum())
    negative_total = float(negative.sum())
    absolute_total = positive_total + negative_total
    interstitial_positive = float(positive[~owned].sum())
    cavity_positive = float(positive[partition.cavity_indices].sum())
    denominator = positive_total if positive_total > 0 else 1.0
    return {
        "signed_integral": float(signed.sum()),
        "positive_integral": positive_total,
        "negative_magnitude_integral": negative_total,
        "absolute_integral": absolute_total,
        "atom_signed_integrals": atom_signed,
        "atom_positive_integrals": atom_positive,
        "atom_negative_magnitude_integrals": atom_negative,
        "molecule_signed_integrals": molecule_signed,
        "molecule_positive_integrals": molecule_positive,
        "molecule_negative_magnitude_integrals": molecule_negative,
        "interstitial_positive_integral": interstitial_positive,
        "interstitial_positive_fraction": interstitial_positive / denominator,
        "cavity_probe_positive_integral": cavity_positive,
        "cavity_probe_positive_fraction": cavity_positive / denominator,
    }


def _periodic_distance(
    first: Sequence[float], second: Sequence[float], cell_angstrom: NDArray[np.float64]
) -> float:
    vector = minimum_image_vectors(
        np.asarray(first, dtype=float).reshape(1, 3),
        np.asarray(second, dtype=float).reshape(1, 3),
        cell_angstrom,
    )[0, 0]
    return float(np.linalg.norm(vector))


def _nearest_real_atom_surface(
    point_angstrom: NDArray[np.float64],
    elements: Sequence[str],
    positions_angstrom: NDArray[np.float64],
    cell_angstrom: NDArray[np.float64],
) -> dict[str, Any]:
    real_indices = [index for index, element in enumerate(elements) if element.upper() != "GH"]
    if not real_indices:
        raise ValueError("localization analysis requires at least one real atom")
    reference = np.asarray(point_angstrom, dtype=float).reshape(1, 3)
    positions = np.asarray(positions_angstrom, dtype=float)[real_indices]
    vectors = minimum_image_vectors(reference, positions, cell_angstrom)[0]
    distances = np.linalg.norm(vectors, axis=1)
    clearances = []
    for local_index, atom_index in enumerate(real_indices):
        element = elements[atom_index].upper()
        if element not in VDW_RADII_ANGSTROM:
            raise ValueError(f"unsupported element for centroid clearance: {element!r}")
        clearances.append(float(distances[local_index] - VDW_RADII_ANGSTROM[element]))
    selected = int(np.argmin(clearances))
    atom_index = real_indices[selected]
    return {
        "atom_index_1based": atom_index + 1,
        "element": elements[atom_index],
        "center_distance_angstrom": float(distances[selected]),
        "vdw_surface_clearance_angstrom": float(clearances[selected]),
    }


def state_localization_record(
    cube: CubeData,
    elements: Sequence[str],
    positions_angstrom: NDArray[np.float64],
    cell_angstrom: NDArray[np.float64],
    topology: MolecularTopology,
    cavity_center_angstrom: NDArray[np.float64],
    settings: Mapping[str, Any],
) -> dict[str, Any]:
    metrics = analyze_spin_density(cube)
    scale = length_scale_to_angstrom(cube)
    centroid = (
        np.asarray(metrics.centroid, dtype=float) * scale
        if metrics.centroid is not None
        else None
    )
    partition = build_geometric_partition(
        cube,
        elements,
        positions_angstrom,
        cavity_center_angstrom,
        vdw_scale=float(settings["vdw_region_scale"]),
        cavity_probe_radius_angstrom=float(settings["cavity_probe_radius_angstrom"]),
    )
    integrated = integrate_partitioned_field(
        cube, cube.values, partition, topology.molecule_ids
    )
    positive_total = float(integrated["positive_integral"])
    denominator = positive_total if positive_total > 0 else 1.0
    li_indices = [index for index, element in enumerate(elements) if element.upper() == "LI"]
    if len(li_indices) > 1:
        raise ValueError(f"candidate may contain at most one Li atom, found {len(li_indices)}")
    li_index = li_indices[0] if li_indices else None
    atom_positive = np.asarray(integrated.pop("atom_positive_integrals"), dtype=float)
    atom_negative = np.asarray(
        integrated.pop("atom_negative_magnitude_integrals"), dtype=float
    )
    atom_signed = np.asarray(integrated.pop("atom_signed_integrals"), dtype=float)
    molecule_positive = np.asarray(integrated.pop("molecule_positive_integrals"), dtype=float)
    molecule_negative = np.asarray(
        integrated.pop("molecule_negative_magnitude_integrals"), dtype=float
    )
    molecule_signed = np.asarray(integrated.pop("molecule_signed_integrals"), dtype=float)
    molecule_records: list[dict[str, Any]] = []
    for molecule in topology.molecules:
        molecule_id = int(molecule["molecule_id"])
        array_index = molecule_id - 1
        molecule_records.append(
            {
                **molecule,
                "signed_spin_electrons": float(molecule_signed[array_index]),
                "positive_spin_electrons": float(molecule_positive[array_index]),
                "negative_spin_magnitude_electrons": float(molecule_negative[array_index]),
                "positive_spin_fraction": float(molecule_positive[array_index] / denominator),
            }
        )
    molecule_records.sort(
        key=lambda record: (-record["positive_spin_electrons"], record["molecule_id"])
    )
    solvent_indices = [
        index for index, molecule_id in enumerate(topology.molecule_ids) if molecule_id >= 0
    ]
    max_atom_index = max(solvent_indices, key=lambda index: atom_positive[index])
    max_molecule_fraction = (
        float(molecule_records[0]["positive_spin_fraction"]) if molecule_records else 0.0
    )
    li_fraction = float(atom_positive[li_index] / denominator) if li_index is not None else None
    centroid_probe_fraction = None
    nearest_atom = None
    if centroid is not None:
        centroid_partition = build_geometric_partition(
            cube,
            elements,
            positions_angstrom,
            centroid,
            vdw_scale=float(settings["vdw_region_scale"]),
            cavity_probe_radius_angstrom=float(settings["cavity_probe_radius_angstrom"]),
        )
        centroid_integrated = integrate_partitioned_field(
            cube, cube.values, centroid_partition, topology.molecule_ids
        )
        centroid_probe_fraction = float(centroid_integrated["cavity_probe_positive_fraction"])
        nearest_atom = _nearest_real_atom_surface(
            centroid, elements, positions_angstrom, cell_angstrom
        )
    proxy_flags = {
        "li_centered": (
            li_fraction is not None
            and li_fraction >= float(settings["li_positive_spin_fraction_threshold"])
        ),
        "single_solvent_molecule": max_molecule_fraction
        >= float(settings["molecular_positive_spin_fraction_threshold"]),
        "ghost_cavity_centered": float(integrated["cavity_probe_positive_fraction"])
        >= float(settings["cavity_positive_spin_fraction_threshold"]),
        "interstitial": float(integrated["interstitial_positive_fraction"])
        >= float(settings["interstitial_positive_spin_fraction_threshold"]),
        "centroid_in_geometric_void": (
            nearest_atom is not None and float(nearest_atom["vdw_surface_clearance_angstrom"]) > 0.0
        ),
    }
    return {
        "spin_density": {
            **metrics.as_dict(),
            "centroid_angstrom": centroid.tolist() if centroid is not None else None,
            "centroid_defined": centroid is not None,
            "radius_angstrom": (
                float(metrics.radius * scale) if metrics.radius is not None else None
            ),
            "negative_magnitude_integral": float(integrated["negative_magnitude_integral"]),
            "absolute_integral": float(integrated["absolute_integral"]),
            "centroid_to_li_angstrom": (
                _periodic_distance(centroid, positions_angstrom[li_index], cell_angstrom)
                if centroid is not None and li_index is not None
                else None
            ),
            "centroid_to_ghost_cavity_angstrom": (
                _periodic_distance(centroid, cavity_center_angstrom, cell_angstrom)
                if centroid is not None
                else None
            ),
            "nearest_real_atom": nearest_atom,
            "centroid_probe_positive_spin_fraction": centroid_probe_fraction,
        },
        "geometric_partition": {
            "definition": "nearest scaled van-der-Waals sphere; outside union is interstitial",
            "vdw_region_scale": float(settings["vdw_region_scale"]),
            "cavity_probe_radius_angstrom": float(settings["cavity_probe_radius_angstrom"]),
            "li_atom_present": li_index is not None,
            "li_signed_spin_electrons": (
                float(atom_signed[li_index]) if li_index is not None else None
            ),
            "li_positive_spin_electrons": (
                float(atom_positive[li_index]) if li_index is not None else None
            ),
            "li_negative_spin_magnitude_electrons": (
                float(atom_negative[li_index]) if li_index is not None else None
            ),
            "li_positive_spin_fraction": li_fraction,
            "maximum_solvent_atom": {
                "atom_index_1based": max_atom_index + 1,
                "element": elements[max_atom_index],
                "signed_spin_electrons": float(atom_signed[max_atom_index]),
                "positive_spin_electrons": float(atom_positive[max_atom_index]),
                "positive_spin_fraction": float(atom_positive[max_atom_index] / denominator),
            },
            "maximum_solvent_molecule_positive_spin_fraction": max_molecule_fraction,
            "interstitial_positive_spin_electrons": float(
                integrated["interstitial_positive_integral"]
            ),
            "interstitial_positive_spin_fraction": float(
                integrated["interstitial_positive_fraction"]
            ),
            "ghost_cavity_probe_positive_spin_electrons": float(
                integrated["cavity_probe_positive_integral"]
            ),
            "ghost_cavity_probe_positive_spin_fraction": float(
                integrated["cavity_probe_positive_fraction"]
            ),
        },
        "localization_proxy_flags": proxy_flags,
        "molecules": molecule_records,
    }


def density_difference_record(
    li0_cube: CubeData,
    separated_cube: CubeData,
    elements: Sequence[str],
    positions_angstrom: NDArray[np.float64],
    topology: MolecularTopology,
    cavity_center_angstrom: NDArray[np.float64],
    settings: Mapping[str, Any],
) -> dict[str, Any]:
    tolerance = float(settings["cube_grid_tolerance_angstrom"])
    problems = cube_compatibility_problems(
        li0_cube, separated_cube, tolerance_angstrom=tolerance
    )
    if problems:
        return {"ready": False, "problems": problems}
    difference = np.asarray(separated_cube.values - li0_cube.values, dtype=float)
    partition = build_geometric_partition(
        li0_cube,
        elements,
        positions_angstrom,
        cavity_center_angstrom,
        vdw_scale=float(settings["vdw_region_scale"]),
        cavity_probe_radius_angstrom=float(settings["cavity_probe_radius_angstrom"]),
    )
    integrated = integrate_partitioned_field(
        li0_cube, difference, partition, topology.molecule_ids
    )
    signed = float(integrated["signed_integral"])
    positive = float(integrated["positive_integral"])
    negative = float(integrated["negative_magnitude_integral"])
    absolute = float(integrated["absolute_integral"])
    conservation_tolerance = float(settings["density_difference_charge_tolerance_electrons"])
    if abs(signed) > conservation_tolerance:
        problems.append(
            f"density-difference integral {signed:.6g} e exceeds {conservation_tolerance:.6g} e"
        )
    molecule_positive = np.asarray(integrated.pop("molecule_positive_integrals"), dtype=float)
    molecule_negative = np.asarray(
        integrated.pop("molecule_negative_magnitude_integrals"), dtype=float
    )
    molecule_signed = np.asarray(integrated.pop("molecule_signed_integrals"), dtype=float)
    atom_positive = np.asarray(integrated.pop("atom_positive_integrals"), dtype=float)
    atom_negative = np.asarray(
        integrated.pop("atom_negative_magnitude_integrals"), dtype=float
    )
    atom_signed = np.asarray(integrated.pop("atom_signed_integrals"), dtype=float)
    li_index = next(index for index, element in enumerate(elements) if element.upper() == "LI")
    molecules: list[dict[str, Any]] = []
    for molecule in topology.molecules:
        array_index = int(molecule["molecule_id"]) - 1
        molecules.append(
            {
                **molecule,
                "density_change_electrons": float(molecule_signed[array_index]),
                "density_accumulation_electrons": float(molecule_positive[array_index]),
                "density_depletion_electrons": float(molecule_negative[array_index]),
            }
        )
    molecules.sort(
        key=lambda record: (-abs(record["density_change_electrons"]), record["molecule_id"])
    )
    return {
        "ready": not problems,
        "problems": problems,
        "definition": "rho(li_plus_e_diabatic) - rho(li0_diabatic)",
        "signed_integral_electrons": signed,
        "accumulation_integral_electrons": positive,
        "depletion_magnitude_integral_electrons": negative,
        "absolute_integral_electrons": absolute,
        "rearranged_electrons_half_l1": 0.5 * absolute,
        "charge_conservation_tolerance_electrons": conservation_tolerance,
        "li_region": {
            "density_change_electrons": float(atom_signed[li_index]),
            "density_accumulation_electrons": float(atom_positive[li_index]),
            "density_depletion_electrons": float(atom_negative[li_index]),
        },
        "interstitial_region": {
            "density_accumulation_electrons": float(
                integrated["interstitial_positive_integral"]
            ),
            "fraction_of_all_accumulation": float(
                integrated["interstitial_positive_fraction"]
            ),
        },
        "ghost_cavity_probe": {
            "radius_angstrom": float(settings["cavity_probe_radius_angstrom"]),
            "density_accumulation_electrons": float(
                integrated["cavity_probe_positive_integral"]
            ),
            "fraction_of_all_accumulation": float(
                integrated["cavity_probe_positive_fraction"]
            ),
        },
        "molecules": molecules,
    }
