from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .localization_analysis import MolecularTopology


def _minimum_image_vector(
    point: NDArray[np.float64],
    reference: NDArray[np.float64],
    cell_angstrom: NDArray[np.float64],
) -> NDArray[np.float64]:
    inverse = np.linalg.inv(cell_angstrom)
    fractional = (np.asarray(point, dtype=float) - np.asarray(reference, dtype=float)) @ inverse
    fractional -= np.rint(fractional)
    return fractional @ cell_angstrom


def periodic_distance(
    first: Sequence[float],
    second: Sequence[float],
    cell_angstrom: NDArray[np.float64],
) -> float:
    vector = _minimum_image_vector(
        np.asarray(first, dtype=float), np.asarray(second, dtype=float), cell_angstrom
    )
    return float(np.linalg.norm(vector))


def molecular_heavy_atom_centers(
    elements: Sequence[str],
    positions_angstrom: NDArray[np.float64],
    cell_angstrom: NDArray[np.float64],
    topology: MolecularTopology,
) -> tuple[dict[str, Any], ...]:
    """Return wrapped heavy-atom geometric centers for inferred solvent molecules."""

    positions = np.asarray(positions_angstrom, dtype=float)
    cell = np.asarray(cell_angstrom, dtype=float)
    inverse = np.linalg.inv(cell)
    records: list[dict[str, Any]] = []
    for molecule in topology.molecules:
        indices = [int(value) - 1 for value in molecule["atom_indices_1based"]]
        heavy = [index for index in indices if str(elements[index]).upper() != "H"]
        selected = heavy or indices
        reference = positions[selected[0]]
        unwrapped = np.vstack(
            [
                reference + _minimum_image_vector(positions[index], reference, cell)
                for index in selected
            ]
        )
        center = np.mean(unwrapped, axis=0)
        fractional = center @ inverse
        center = (fractional - np.floor(fractional)) @ cell
        records.append(
            {
                "molecule_id": int(molecule["molecule_id"]),
                "component": molecule.get("component"),
                "heavy_atom_count": len(heavy),
                "center_angstrom": center.tolist(),
            }
        )
    return tuple(records)


def local_composition(
    point_angstrom: Sequence[float],
    molecule_centers: Sequence[Mapping[str, Any]],
    cell_angstrom: NDArray[np.float64],
    radii_angstrom: Sequence[float],
    component_counts: Mapping[str, int],
    *,
    amine_component: str | None,
) -> dict[str, Any]:
    """Count molecule centers in periodic shells and report amine enrichment."""

    radii = sorted({float(value) for value in radii_angstrom})
    if not radii or radii[0] <= 0:
        raise ValueError("local-composition shell radii must be positive")
    total_bulk = sum(int(value) for value in component_counts.values())
    if total_bulk <= 0:
        raise ValueError("component counts must contain solvent molecules")
    amine_bulk = int(component_counts.get(amine_component, 0)) if amine_component else 0
    bulk_fraction = amine_bulk / total_bulk
    point = np.asarray(point_angstrom, dtype=float)
    distances = [
        (
            periodic_distance(point, record["center_angstrom"], cell_angstrom),
            str(record.get("component")),
        )
        for record in molecule_centers
    ]
    shells: list[dict[str, Any]] = []
    for radius in radii:
        selected = [component for distance, component in distances if distance <= radius]
        counts = {
            component: selected.count(component)
            for component in sorted({str(value) for value in component_counts})
        }
        total = len(selected)
        amine_count = selected.count(str(amine_component)) if amine_component else 0
        fraction = amine_count / total if total else None
        shells.append(
            {
                "radius_angstrom": radius,
                "molecule_count": total,
                "component_counts": counts,
                "amine_count": amine_count,
                "amine_mole_fraction": fraction,
                "bulk_amine_mole_fraction": bulk_fraction,
                "amine_fraction_excess": fraction - bulk_fraction if fraction is not None else None,
                "amine_enrichment_ratio": (
                    fraction / bulk_fraction
                    if fraction is not None and bulk_fraction > 0
                    else None
                ),
            }
        )
    return {
        "point_angstrom": point.tolist(),
        "amine_component": amine_component,
        "bulk_component_counts": {
            str(name): int(value) for name, value in sorted(component_counts.items())
        },
        "bulk_amine_mole_fraction": bulk_fraction,
        "shells": shells,
    }


def select_preferential_sites(
    sites: Sequence[Mapping[str, Any]],
    molecule_centers: Sequence[Mapping[str, Any]],
    cell_angstrom: NDArray[np.float64],
    component_counts: Mapping[str, int],
    *,
    amine_component: str | None,
    shell_radii_angstrom: Sequence[float],
    selection_radius_angstrom: float,
    minimum_separation_angstrom: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Choose deterministic amine-rich and amine-poor void-basis seeds."""

    if len(sites) < 2:
        raise ValueError("at least two void sites are required")
    if selection_radius_angstrom <= 0 or minimum_separation_angstrom <= 0:
        raise ValueError("selection radius and seed separation must be positive")
    radii = sorted({float(value) for value in shell_radii_angstrom})
    if selection_radius_angstrom not in radii:
        raise ValueError("selection radius must be one of the reported shell radii")
    scored: list[dict[str, Any]] = []
    for site in sites:
        point = np.asarray(site["cartesian_angstrom"], dtype=float)
        composition = local_composition(
            point,
            molecule_centers,
            cell_angstrom,
            radii,
            component_counts,
            amine_component=amine_component,
        )
        selection_shell = next(
            shell
            for shell in composition["shells"]
            if float(shell["radius_angstrom"]) == float(selection_radius_angstrom)
        )
        scored.append(
            {
                "site": dict(site),
                "local_composition": composition,
                "selection_amine_mole_fraction": selection_shell["amine_mole_fraction"],
                "selection_molecule_count": int(selection_shell["molecule_count"]),
            }
        )
    scored = [record for record in scored if int(record["selection_molecule_count"]) > 0]
    if len(scored) < 2:
        raise ValueError("fewer than two void sites have solvent molecules in the selection shell")
    scored.sort(key=lambda record: int(record["site"]["rank"]))

    def score(record: Mapping[str, Any]) -> float:
        value = record["selection_amine_mole_fraction"]
        return float(value) if value is not None else -1.0

    rich = min(
        scored,
        key=lambda record: (
            -score(record),
            -int(record["selection_molecule_count"]),
            int(record["site"]["rank"]),
        ),
    )
    eligible = [
        record
        for record in scored
        if record is not rich
        and periodic_distance(
            rich["site"]["cartesian_angstrom"],
            record["site"]["cartesian_angstrom"],
            cell_angstrom,
        )
        >= minimum_separation_angstrom
    ]
    if not eligible:
        raise ValueError("no second void site satisfies the minimum seed separation")
    poor = min(
        eligible,
        key=lambda record: (
            score(record),
            int(record["selection_molecule_count"]),
            int(record["site"]["rank"]),
        ),
    )
    fractions = {score(record) for record in scored}
    degenerate = len(fractions) == 1

    def finalize(role: str, record: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "seed_role": role,
            "site": record["site"],
            "local_composition": record["local_composition"],
            "selection_radius_angstrom": float(selection_radius_angstrom),
            "selection_amine_mole_fraction": record["selection_amine_mole_fraction"],
            "composition_degenerate_reference": degenerate,
        }

    return finalize("eda_rich", rich), finalize("eda_poor", poor)
