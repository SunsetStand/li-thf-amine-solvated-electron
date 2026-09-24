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
    target_electrons: float | None = None,
) -> None:
    """Render the deliberately low-cost Stage-B numerical cDFT smoke input."""

    if li_atom_index <= 0:
        raise ValueError("CP2K atom indices are one-based")
    if method.get("scientific_status") != "NUMERICAL_SMOKE_ONLY":
        raise ValueError("Stage-B smoke method must remain explicitly non-production")
    li_valence = int(method["li_pseudopotential_valence_electrons"])
    configured_target = float(method["li_target_valence_electrons"])
    if li_valence - configured_target != 1.0:
        raise ValueError("Li+ cDFT requires target = pseudopotential valence - 1")
    if target_electrons is None:
        target = configured_target
    else:
        target = float(target_electrons)
        allowed_targets = {float(li_valence), float(li_valence - 1)}
        if target not in allowed_targets:
            raise ValueError(
                "Stage-B mechanism smoke target must represent Li0 or Li+ "
                f"({sorted(allowed_targets)} valence electrons)"
            )
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


def render_stage_b2_intrinsic_cp2k(
    template_path: str | Path,
    output_path: str | Path,
    *,
    project: str,
    coordinates_path: str | Path,
    cell_path: str | Path,
    method: Mapping[str, Any],
) -> None:
    """Render a Li-free, unconstrained excess-electron numerical smoke input."""

    if method.get("scientific_status") != "NUMERICAL_VERTICAL_ELECTRON_ONLY_SMOKE":
        raise ValueError("Stage-B2 intrinsic method must remain explicitly smoke-only")
    if int(method.get("charge", 0)) != -1 or int(method.get("multiplicity", 0)) != 2:
        raise ValueError("Stage-B2 intrinsic smoke requires a charge -1 doublet")
    substitutions = {
        "project": project,
        "coordinates_path": Path(coordinates_path).resolve().as_posix(),
        "cell_path": Path(cell_path).resolve().as_posix(),
        "charge": int(method["charge"]),
        "multiplicity": int(method["multiplicity"]),
        "basis_set": method["basis_set"],
        "ghost_basis_set": method["ghost_basis_set"],
        "potential": method["potential"],
        "cutoff_ry": method["cutoff_ry"],
        "rel_cutoff_ry": method["rel_cutoff_ry"],
        "eps_scf": method["eps_scf"],
        "max_scf": int(method["max_scf"]),
        "cube_stride": int(method["cube_stride"]),
    }
    rendered = _load_template(template_path).substitute(substitutions)
    upper = rendered.upper()
    if "&CDFT" in upper or "&CONSTRAINT" in upper:
        raise ValueError("Stage-B2 intrinsic smoke must not contain localization constraints")
    if "&KIND LI" in upper:
        raise ValueError("Stage-B2 intrinsic smoke must not define a Li kind")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered.rstrip() + "\n", encoding="utf-8")


def render_stage_b2c_preferential_cp2k(
    template_path: str | Path,
    output_path: str | Path,
    *,
    project: str,
    coordinates_path: str | Path,
    cell_path: str | Path,
    method: Mapping[str, Any],
    state: Mapping[str, Any],
) -> None:
    """Render one member of a fixed-geometry neutral/anion Stage-B2C pair."""

    if method.get("scientific_status") != "NUMERICAL_PREFERENTIAL_SOLVATION_SMOKE_ONLY":
        raise ValueError("Stage-B2C method must remain explicitly smoke-only")
    state_id = str(state.get("id", ""))
    allowed = {
        "neutral": (0, 1, False),
        "anion": (-1, 2, True),
    }
    charge = int(state.get("charge", 99))
    multiplicity = int(state.get("multiplicity", 99))
    uks = bool(state.get("uks"))
    if state_id not in allowed or (charge, multiplicity, uks) != allowed[state_id]:
        raise ValueError("Stage-B2C states must be neutral singlet and charge -1 doublet")
    density_print = ""
    if state_id == "anion":
        stride = int(method["cube_stride"])
        density_print = f"""      &E_DENSITY_CUBE
        STRIDE {stride} {stride} {stride}
      &END E_DENSITY_CUBE"""
    substitutions = {
        "project": project,
        "coordinates_path": Path(coordinates_path).resolve().as_posix(),
        "cell_path": Path(cell_path).resolve().as_posix(),
        "charge": charge,
        "multiplicity": multiplicity,
        "uks": "TRUE" if uks else "FALSE",
        "basis_set": method["basis_set"],
        "ghost_basis_set": method["ghost_basis_set"],
        "potential": method["potential"],
        "cutoff_ry": method["cutoff_ry"],
        "rel_cutoff_ry": method["rel_cutoff_ry"],
        "eps_scf": method["eps_scf"],
        "max_scf": int(method["max_scf"]),
        "density_print_block": density_print,
    }
    rendered = _load_template(template_path).substitute(substitutions)
    upper = rendered.upper()
    if "&CDFT" in upper or "&CONSTRAINT" in upper:
        raise ValueError("Stage-B2C smoke must not contain localization constraints")
    if "&KIND LI" in upper:
        raise ValueError("Stage-B2C smoke must not define a Li kind")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered.rstrip() + "\n", encoding="utf-8")


def render_stage_b3_amine_series_cp2k(
    template_path: str | Path,
    output_path: str | Path,
    *,
    project: str,
    coordinates_path: str | Path,
    cell_path: str | Path,
    method: Mapping[str, Any],
    state: Mapping[str, Any],
) -> None:
    """Render a fixed-nuclei Stage-B3 neutral/anion screening state."""

    if method.get("scientific_status") != (
        "NUMERICAL_AMINE_SERIES_NH_ASSOCIATION_SCREEN_ONLY"
    ):
        raise ValueError("Stage-B3 method must remain explicitly screening-only")
    compatible = dict(method)
    compatible["scientific_status"] = "NUMERICAL_PREFERENTIAL_SOLVATION_SMOKE_ONLY"
    render_stage_b2c_preferential_cp2k(
        template_path,
        output_path,
        project=project,
        coordinates_path=coordinates_path,
        cell_path=cell_path,
        method=compatible,
        state=state,
    )


def render_stage_b2c_method_benchmark_cp2k(
    template_path: str | Path,
    output_path: str | Path,
    *,
    project: str,
    coordinates_path: str | Path,
    cell_path: str | Path,
    method: Mapping[str, Any],
    state: Mapping[str, Any],
) -> None:
    """Render one fixed-nuclei Stage-B2C method-benchmark state."""

    if method.get("scientific_status") != (
        "NUMERICAL_METHOD_SENSITIVITY_BENCHMARK_ONLY"
    ):
        raise ValueError("Stage-B2C benchmark method must remain explicitly diagnostic")
    state_id = str(state.get("id", ""))
    allowed = {
        "neutral": (0, 1, False),
        "anion": (-1, 2, True),
    }
    charge = int(state.get("charge", 99))
    multiplicity = int(state.get("multiplicity", 99))
    uks = bool(state.get("uks"))
    if state_id not in allowed or (charge, multiplicity, uks) != allowed[state_id]:
        raise ValueError("Stage-B2C benchmark states must be neutral singlet and anion doublet")

    family = str(method.get("xc_family", "")).upper()
    admm = bool(method.get("admm"))
    fraction = float(method.get("exact_exchange_fraction", 0.0))
    basis_files = "    BASIS_SET_FILE_NAME BASIS_MOLOPT"
    aux_fit_line = ""
    admm_block = ""
    if family == "PBE":
        if admm or fraction != 0.0:
            raise ValueError("PBE benchmark variant must not enable exact exchange or ADMM")
        xc_lines = [
            "      &XC_FUNCTIONAL PBE",
            "      &END XC_FUNCTIONAL",
        ]
        dispersion_reference = "PBE"
    elif family == "PBE0":
        if not admm or not 0.0 < fraction < 1.0:
            raise ValueError("PBE0 benchmark variant requires ADMM and partial exact exchange")
        basis_files += "\n    BASIS_SET_FILE_NAME BASIS_ADMM_MOLOPT"
        aux_basis = str(method["aux_basis_set"])
        aux_fit_line = f"      BASIS_SET AUX_FIT {aux_basis}"
        admm_block = """    &AUXILIARY_DENSITY_MATRIX_METHOD
      ADMM_TYPE ADMMS
      EXCH_CORRECTION_FUNC PBEX
    &END AUXILIARY_DENSITY_MATRIX_METHOD"""
        xc_lines = [
            "      &XC_FUNCTIONAL PBE",
            "      &END XC_FUNCTIONAL",
            "      &HF",
            f"        FRACTION {fraction}",
            "        &SCREENING",
            f"          EPS_SCHWARZ {method['eps_schwarz']}",
            "          SCREEN_ON_INITIAL_P TRUE",
            "        &END SCREENING",
            "        &INTERACTION_POTENTIAL",
            "          POTENTIAL_TYPE TRUNCATED",
            f"          CUTOFF_RADIUS {method['hfx_cutoff_angstrom']}",
            "          T_C_G_DATA t_c_g.dat",
            "        &END INTERACTION_POTENTIAL",
            "        &MEMORY",
            "          MAX_MEMORY 3000",
            "          EPS_STORAGE_SCALING 0.1",
            "        &END MEMORY",
            "      &END HF",
        ]
        dispersion_reference = "PBE0"
    else:
        raise ValueError(f"unsupported Stage-B2C benchmark XC family: {family!r}")
    xc_lines.extend(
        [
            "      &VDW_POTENTIAL",
            "        POTENTIAL_TYPE PAIR_POTENTIAL",
            "        &PAIR_POTENTIAL",
            "          TYPE DFTD3(BJ)",
            "          PARAMETER_FILE_NAME dftd3.dat",
            f"          REFERENCE_FUNCTIONAL {dispersion_reference}",
            "          R_CUTOFF 15.0",
            "        &END PAIR_POTENTIAL",
            "      &END VDW_POTENTIAL",
        ]
    )

    density_print = ""
    if state_id == "anion":
        stride = int(method["cube_stride"])
        density_print = f"""      &E_DENSITY_CUBE
        STRIDE {stride} {stride} {stride}
      &END E_DENSITY_CUBE"""
    substitutions = {
        "project": project,
        "coordinates_path": Path(coordinates_path).resolve().as_posix(),
        "cell_path": Path(cell_path).resolve().as_posix(),
        "charge": charge,
        "multiplicity": multiplicity,
        "uks": "TRUE" if uks else "FALSE",
        "basis_file_block": basis_files,
        "basis_set": method["basis_set"],
        "ghost_basis_set": method["ghost_basis_set"],
        "aux_fit_line": aux_fit_line,
        "potential": method["potential"],
        "cutoff_ry": method["cutoff_ry"],
        "rel_cutoff_ry": method["rel_cutoff_ry"],
        "eps_scf": method["eps_scf"],
        "max_scf": int(method["max_scf"]),
        "admm_block": admm_block,
        "xc_block": "\n".join(xc_lines),
        "density_print_block": density_print,
    }
    rendered = _load_template(template_path).substitute(substitutions)
    upper = rendered.upper()
    if "&CDFT" in upper or "&CONSTRAINT" in upper or "&KIND LI" in upper:
        raise ValueError("Stage-B2C benchmark must not contain Li or localization constraints")
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
    if box_angstrom <= 0 or thf_count < 0 or amine_count < 0:
        raise ValueError("box must be positive and molecule counts non-negative")
    if thf_count + amine_count <= 0:
        raise ValueError("at least one molecule is required")
    lines = [
        "tolerance 2.0",
        "filetype pdb",
        f"output {output_pdb}",
        f"seed {seed}",
        "add_box_sides 1.0",
    ]
    if thf_count:
        lines.extend(
            [
                "",
                f"structure {thf_structure}",
                f"  number {thf_count}",
                "  inside box 0.0 0.0 0.0 "
                f"{box_angstrom:.6f} {box_angstrom:.6f} {box_angstrom:.6f}",
                "end structure",
            ]
        )
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
