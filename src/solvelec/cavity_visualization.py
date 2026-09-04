from __future__ import annotations

import html
from pathlib import Path
from typing import Any

import numpy as np

from .trajectory import minimum_image_vectors

ELEMENT_COLORS = {
    "H": "#f8fafc",
    "C": "#64748b",
    "N": "#2563eb",
    "O": "#dc2626",
    "LI": "#7c3aed",
}
ELEMENT_RADII = {"H": 4.0, "C": 7.0, "N": 8.0, "O": 8.0, "LI": 10.0}


def cavity_centered_positions(
    positions: np.ndarray,
    residue_ids: np.ndarray,
    heavy_indices: np.ndarray,
    cavity_center: np.ndarray,
    cell: np.ndarray,
) -> np.ndarray:
    """Translate intact residues to their nearest image around a cavity center."""

    coordinates = np.asarray(positions, dtype=float)
    residues = np.asarray(residue_ids)
    heavy = set(int(value) for value in np.asarray(heavy_indices, dtype=int))
    center = np.asarray(cavity_center, dtype=float)
    localized = np.array(coordinates, copy=True)
    for residue_id in np.unique(residues):
        indices = np.flatnonzero(residues == residue_id)
        anchor_indices = [index for index in indices if int(index) in heavy]
        if not anchor_indices:
            anchor_indices = list(indices)
        anchor_position = coordinates[[anchor_indices[0]]]
        # Reassemble the residue before taking its centroid. A plain mean is
        # wrong when a wrapped molecule straddles a periodic boundary.
        residue_positions = anchor_position + minimum_image_vectors(
            anchor_position, coordinates[indices], cell
        )[0]
        heavy_offsets = [
            offset for offset, atom_index in enumerate(indices) if int(atom_index) in heavy
        ]
        if not heavy_offsets:
            heavy_offsets = list(range(len(indices)))
        residue_center = np.mean(residue_positions[heavy_offsets], axis=0, keepdims=True)
        nearest_vector = minimum_image_vectors(center[None, :], residue_center, cell)[0, 0]
        localized[indices] = residue_positions - residue_center + center + nearest_vector
    return localized - center


def local_residue_mask(
    centered_positions: np.ndarray,
    residue_ids: np.ndarray,
    heavy_indices: np.ndarray,
    radius_angstrom: float,
) -> np.ndarray:
    if radius_angstrom <= 0:
        raise ValueError("visualization radius must be positive")
    residues = np.asarray(residue_ids)
    positions = np.asarray(centered_positions, dtype=float)
    heavy = set(int(value) for value in np.asarray(heavy_indices, dtype=int))
    selected_residues: set[int] = set()
    for index in heavy:
        if float(np.linalg.norm(positions[index])) <= radius_angstrom:
            selected_residues.add(int(residues[index]))
    return np.asarray([int(value) in selected_residues for value in residues], dtype=bool)


def write_local_pdb(path: Path, atoms: list[dict[str, Any]]) -> dict[int, int]:
    """Write a cavity-centered local cluster and return atom-index to PDB-serial mapping."""

    path.parent.mkdir(parents=True, exist_ok=True)
    serials: dict[int, int] = {}
    lines = ["REMARK cavity-centered solvent shell; cavity center is at 0 0 0"]
    for serial, atom in enumerate(atoms, start=1):
        serials[int(atom["index"])] = serial
        x, y, z = (float(value) for value in atom["position_angstrom"])
        atom_name = str(atom["name"])[:4]
        residue_name = str(atom["resname"])[:3]
        residue_id = int(atom["resid"]) % 10000
        element = str(atom["element"]).upper()[:2]
        lines.append(
            f"HETATM{serial:5d} {atom_name:<4s} {residue_name:>3s} A{residue_id:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
        )
    lines.extend(["TER", "END"])
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(path)
    return serials


def write_pymol_script(
    path: Path,
    pdb_path: Path,
    cavity_radius_angstrom: float,
    bonds: list[dict[str, Any]],
    serials: dict[int, int],
) -> None:
    lines = [
        f"load {pdb_path.name}, shell",
        "hide everything, all",
        "show sticks, shell",
        "set stick_radius, 0.12, shell",
        "color gray70, elem C",
        "color red, elem O",
        "color marine, elem N",
        "color white, elem H",
        "pseudoatom cavity_center, pos=[0.0, 0.0, 0.0]",
        "show spheres, cavity_center",
        f"set sphere_scale, {cavity_radius_angstrom:.4f}, cavity_center",
        "set sphere_transparency, 0.65, cavity_center",
        "color cyan, cavity_center",
    ]
    for number, bond in enumerate(bonds, start=1):
        hydrogen = serials.get(int(bond["hydrogen_index"]))
        acceptor = serials.get(int(bond["acceptor_index"]))
        if hydrogen is None or acceptor is None:
            continue
        name = f"hbond_{number:03d}"
        lines.append(f"distance {name}, shell and id {hydrogen}, shell and id {acceptor}")
        color = "orange" if bond["bond_type"] == "eda_nh_eda_n" else "teal"
        lines.append(f"color {color}, {name}")
        lines.append(f"set dash_width, {4 if bond['cavity_bridging'] else 2}, {name}")
    lines.extend(
        [
            "set dash_gap, 0.18",
            "set dash_length, 0.25",
            "bg_color white",
            "orient shell",
            "zoom shell, 2.0",
        ]
    )
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(path)


def _svg_line(x1: float, y1: float, x2: float, y2: float, **attrs: str) -> str:
    values = " ".join(
        f'{name.replace("_", "-")}="{html.escape(value)}"' for name, value in attrs.items()
    )
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" {values}/>'


def render_cavity_hbond_svg(
    path: Path,
    *,
    title: str,
    atoms: list[dict[str, Any]],
    covalent_bonds: list[tuple[int, int]],
    hydrogen_bonds: list[dict[str, Any]],
    cavity_radius_angstrom: float,
    view_radius_angstrom: float,
) -> None:
    """Render three deterministic projections of a cavity and its H-bond shell."""

    width, height = 1800, 720
    panel_width = 560
    panel_height = 560
    top = 105
    gap = 25
    projections = ((0, 1, "xy"), (0, 2, "xz"), (1, 2, "yz"))
    atom_by_index = {int(atom["index"]): atom for atom in atoms}
    scale = panel_height * 0.43 / view_radius_angstrom
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        f'<text x="40" y="45" font-family="Arial,sans-serif" font-size="28" font-weight="700" fill="#0f172a">{html.escape(title)}</text>',
        '<text x="40" y="76" font-family="Arial,sans-serif" font-size="17" fill="#475569">Translucent circle: geometric void proxy · dashed: N-H...N/O · thick dashed: geometrically bridging the void</text>',
    ]
    for panel_index, (axis_x, axis_y, label) in enumerate(projections):
        left = 35 + panel_index * (panel_width + gap)
        center_x = left + panel_width / 2
        center_y = top + panel_height / 2
        lines.append(
            f'<rect x="{left}" y="{top}" width="{panel_width}" height="{panel_height}" rx="18" fill="#ffffff" stroke="#cbd5e1" stroke-width="2"/>'
        )
        lines.append(
            f'<circle cx="{center_x:.2f}" cy="{center_y:.2f}" r="{cavity_radius_angstrom * scale:.2f}" fill="#67e8f9" fill-opacity="0.25" stroke="#0891b2" stroke-width="3"/>'
        )

        def project(position: list[float]) -> tuple[float, float]:
            return (
                center_x + float(position[axis_x]) * scale,
                center_y - float(position[axis_y]) * scale,
            )

        for first, second in covalent_bonds:
            if first not in atom_by_index or second not in atom_by_index:
                continue
            x1, y1 = project(atom_by_index[first]["position_angstrom"])
            x2, y2 = project(atom_by_index[second]["position_angstrom"])
            lines.append(_svg_line(x1, y1, x2, y2, stroke="#94a3b8", stroke_width="3"))
        for bond in hydrogen_bonds:
            hydrogen = atom_by_index.get(int(bond["hydrogen_index"]))
            acceptor = atom_by_index.get(int(bond["acceptor_index"]))
            if hydrogen is None or acceptor is None:
                continue
            x1, y1 = project(hydrogen["position_angstrom"])
            x2, y2 = project(acceptor["position_angstrom"])
            color = "#f97316" if bond["bond_type"] == "eda_nh_eda_n" else "#0f766e"
            width_value = "6" if bond["cavity_bridging"] else "3"
            opacity = "1.0" if bond["cavity_associated"] else "0.28"
            lines.append(
                _svg_line(
                    x1,
                    y1,
                    x2,
                    y2,
                    stroke=color,
                    stroke_width=width_value,
                    stroke_dasharray="12 8",
                    stroke_opacity=opacity,
                )
            )
        depth_axis = ({0, 1, 2} - {axis_x, axis_y}).pop()
        for atom in sorted(atoms, key=lambda value: float(value["position_angstrom"][depth_axis])):
            x, y = project(atom["position_angstrom"])
            element = str(atom["element"]).upper()
            lines.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{ELEMENT_RADII.get(element, 6.0):.1f}" fill="{ELEMENT_COLORS.get(element, "#a3a3a3")}" stroke="#334155" stroke-width="1"/>'
            )
        lines.append(
            f'<text x="{left + 22}" y="{top + 38}" font-family="Arial,sans-serif" font-size="24" font-weight="700" fill="#0f172a">{label.upper()} projection</text>'
        )
    lines.extend(
        [
            '<line x1="50" y1="690" x2="105" y2="690" stroke="#f97316" stroke-width="4" stroke-dasharray="12 8"/>',
            '<text x="115" y="696" font-family="Arial,sans-serif" font-size="18" fill="#334155">EDA N-H...N</text>',
            '<line x1="320" y1="690" x2="375" y2="690" stroke="#0f766e" stroke-width="4" stroke-dasharray="12 8"/>',
            '<text x="385" y="696" font-family="Arial,sans-serif" font-size="18" fill="#334155">EDA N-H...O(THF)</text>',
            '<text x="1660" y="696" text-anchor="end" font-family="Arial,sans-serif" font-size="16" fill="#64748b">Geometric visualization; not electron density</text>',
            "</svg>",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(path)
