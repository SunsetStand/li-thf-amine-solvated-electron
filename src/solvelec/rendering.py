from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from string import Template
from typing import Any

from .state import get_state


def _load_template(path: str | Path) -> Template:
    return Template(Path(path).read_text(encoding="utf-8"))


def render_cp2k(
    template_path: str | Path,
    output_path: str | Path,
    state_name: str,
    project: str,
    coordinates_include: str,
    cell_include: str,
    method: Mapping[str, Any],
    li_atom_index: int,
    constrained: bool,
) -> None:
    state = get_state(state_name)
    cdft = ""
    if constrained:
        if state_name != "solvated_electron":
            raise ValueError("Li+ cDFT constraint is only valid for the solvated-electron state")
        li_valence = int(method["li_pseudopotential_valence_electrons"])
        li_target = float(method["li_target_valence_electrons"])
        if li_valence - li_target != 1.0:
            raise ValueError("Li+ cDFT requires target = pseudopotential valence - 1")
        cdft = f"""      &CDFT
        TYPE_OF_CONSTRAINT BECKE
        ATOMIC_CHARGES TRUE
        STRENGTH 0.0
        TARGET {li_target:.1f}
        &ATOM_GROUP
          ATOMS {li_atom_index}
          COEFF 1.0
          CONSTRAINT_TYPE CHARGE
        &END ATOM_GROUP
        &OUTER_SCF ON
          TYPE CDFT_CONSTRAINT
          EXTRAPOLATION_ORDER 2
          MAX_SCF 20
          EPS_SCF 1.0E-3
          OPTIMIZER NEWTON_LS
          STEP_SIZE -1.0
          &CDFT_OPT ON
            MAX_LS 5
            CONTINUE_LS
            FACTOR_LS 0.5
            JACOBIAN_STEP 1.0E-2
            JACOBIAN_FREQ 1 1
            JACOBIAN_TYPE FD1
            JACOBIAN_RESTART FALSE
          &END CDFT_OPT
        &END OUTER_SCF
        &BECKE_CONSTRAINT
          CUTOFF_TYPE GLOBAL
          GLOBAL_CUTOFF 6.0
          CAVITY_CONFINE TRUE
          CAVITY_SHAPE VDW
          EPS_CAVITY 1.0E-7
          SHOULD_SKIP TRUE
        &END BECKE_CONSTRAINT
      &END CDFT"""
    substitutions = {
        "project": project,
        "charge": state.charge,
        "multiplicity": state.multiplicity,
        "cutoff_ry": method["cutoff_ry"],
        "rel_cutoff_ry": method["rel_cutoff_ry"],
        "basis_set": method["basis_set"],
        "aux_basis_set": method["aux_basis_set"],
        "potential": method["potential"],
        "li_potential": method["li_potential"],
        "exact_exchange_fraction": method["exact_exchange_fraction"],
        "hfx_cutoff_angstrom": method["hfx_cutoff_angstrom"],
        "coordinates_include": coordinates_include,
        "cell_include": cell_include,
        "cdft_block": cdft,
    }
    rendered = _load_template(template_path).substitute(substitutions)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered.rstrip() + "\n", encoding="utf-8")


def render_stage_b_cp2k(
    template_path: str | Path,
    output_path: str | Path,
    *,
    project: str,
    coordinates_path: str | Path,
    cell_path: str | Path,
    method: Mapping[str, Any],
    li_atom_index: int,
) -> None:
    """Render the deliberately low-cost Stage-B numerical cDFT smoke input."""

    if li_atom_index <= 0:
        raise ValueError("CP2K atom indices are one-based")
    if method.get("scientific_status") != "NUMERICAL_SMOKE_ONLY":
        raise ValueError("Stage-B smoke method must remain explicitly non-production")
    li_valence = int(method["li_pseudopotential_valence_electrons"])
    target = float(method["li_target_valence_electrons"])
    if li_valence - target != 1.0:
        raise ValueError("Li+ cDFT requires target = pseudopotential valence - 1")
    substitutions = {
        "project": project,
        "coordinates_path": Path(coordinates_path).resolve().as_posix(),
        "cell_path": Path(cell_path).resolve().as_posix(),
        "basis_set": method["basis_set"],
        "ghost_basis_set": method["ghost_basis_set"],
        "potential": method["potential"],
        "li_potential": method["li_potential"],
        "cutoff_ry": method["cutoff_ry"],
        "rel_cutoff_ry": method["rel_cutoff_ry"],
        "eps_scf": method["eps_scf"],
        "cdft_eps_scf": method["cdft_eps_scf"],
        "max_scf": method["max_scf"],
        "li_atom_index": li_atom_index,
        "li_target_valence_electrons": f"{target:.1f}",
        "cube_stride": int(method["cube_stride"]),
    }
    rendered = _load_template(template_path).substitute(substitutions)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered.rstrip() + "\n", encoding="utf-8")


def render_orca(
    template_path: str | Path,
    output_path: str | Path,
    state_name: str,
    coordinates_xyz: str,
    method: Mapping[str, Any],
) -> None:
    state = get_state(state_name)
    substitutions = {
        "functional": method["functional"],
        "basis_set": method["basis_set"],
        "grid": method["grid"],
        "scf": method["scf"],
        "charge": state.charge,
        "multiplicity": state.multiplicity,
        "coordinates_xyz": coordinates_xyz.rstrip(),
    }
    rendered = _load_template(template_path).substitute(substitutions)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered.rstrip() + "\n", encoding="utf-8")


def render_packmol(
    output_path: str | Path,
    output_pdb: str,
    box_angstrom: float,
    thf_structure: str,
    thf_count: int,
    amine_structure: str | None = None,
    amine_count: int = 0,
    seed: int = 1,
) -> None:
    if box_angstrom <= 0 or thf_count <= 0 or amine_count < 0:
        raise ValueError("box/thf counts must be positive and amine count non-negative")
    lines = [
        "tolerance 2.0",
        "filetype pdb",
        f"output {output_pdb}",
        f"seed {seed}",
        "add_box_sides 1.0",
        "",
        f"structure {thf_structure}",
        f"  number {thf_count}",
        f"  inside box 0.0 0.0 0.0 {box_angstrom:.6f} {box_angstrom:.6f} {box_angstrom:.6f}",
        "end structure",
    ]
    if amine_count:
        if not amine_structure:
            raise ValueError("amine_structure is required when amine_count is non-zero")
        lines.extend(
            [
                "",
                f"structure {amine_structure}",
                f"  number {amine_count}",
                "  inside box 0.0 0.0 0.0 "
                f"{box_angstrom:.6f} {box_angstrom:.6f} {box_angstrom:.6f}",
                "end structure",
            ]
        )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_gromacs_mdp(
    template_path: str | Path,
    output_path: str | Path,
    temperature_k: float,
    pressure_bar: float,
    seed: int,
    *,
    timestep_fs: float = 2.0,
    nsteps: int | None = None,
    trajectory_stride_steps: int = 5000,
) -> None:
    if temperature_k <= 0 or pressure_bar <= 0 or timestep_fs <= 0:
        raise ValueError("temperature, pressure, and timestep must be positive")
    if seed <= 0 or trajectory_stride_steps <= 0:
        raise ValueError("GROMACS seed and trajectory stride must be positive")
    template = _load_template(template_path)
    if "$nsteps" in template.template and (nsteps is None or nsteps <= 0):
        raise ValueError("a positive nsteps value is required by this GROMACS template")
    rendered = template.substitute(
        temperature_k=f"{temperature_k:.2f}",
        pressure_bar=f"{pressure_bar:.6f}",
        seed=seed,
        timestep_ps=f"{timestep_fs / 1000.0:.6f}",
        nsteps=nsteps,
        trajectory_stride_steps=trajectory_stride_steps,
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered.rstrip() + "\n", encoding="utf-8")


def load_spec(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))
