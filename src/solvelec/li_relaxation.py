"""Charge-conserving Li electron transfer with fixed-cell structural relaxation."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from string import Template
from typing import Any

import numpy as np

from .candidates import read_xyz
from .parsers import HARTREE_TO_EV, evaluate_cp2k_cdft_constraint, parse_cp2k_text
from .trajectory import minimum_image_vectors

STATES = ("a", "b")
PAIRS = ("aa", "ba", "bb", "ab")  # electronic state first, nuclear geometry second
STATUS = "FIXED_CELL_RELAXED_LI_ELECTRON_TRANSFER"
VALENCE = {"H": 1, "C": 4, "N": 5, "O": 6, "LI": 3}


def read_json(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value


def write_json(path: str | Path, value: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(output)


def settings_from(methods: dict[str, Any]) -> dict[str, Any]:
    settings = methods["li_relaxation"]
    if settings["scientific_status"] != STATUS:
        raise ValueError("invalid Li relaxation scientific status")
    if settings["systems"] != ["eda_3m", "tmeda_3m"]:
        raise ValueError("Li relaxation is configured for EDA/TMEDA at 3 M")
    if set(settings["source_campaigns"]) != set(settings["systems"]):
        raise ValueError("source campaigns must cover the selected systems")
    if set(settings["source_replicas"]) != set(settings["systems"]) or any(
        not isinstance(value, int) or value < 1 for value in settings["source_replicas"].values()
    ):
        raise ValueError("positive source replica indices must cover the selected systems")
    if set(settings["state_targets"]) != set(STATES):
        raise ValueError("both a and b Li diabatic targets are required")
    targets = [float(settings["state_targets"][state]) for state in STATES]
    if not all(math.isfinite(x) for x in targets) or not 0 < targets[1] < targets[0] <= 3:
        raise ValueError("Li valence population targets must satisfy 0 < b < a <= 3")
    if int(settings["snapshots_per_system"]) < 1:
        raise ValueError("at least one snapshot is required")
    if not 1 <= int(settings["cube_stride"]) <= 4:
        raise ValueError("cube stride must be in [1, 4]")
    if not math.isfinite(float(settings["energy_tolerance_ev"])) or (
        float(settings["energy_tolerance_ev"]) <= 0
    ):
        raise ValueError("energy tolerance must be finite and positive")
    for key in (
        "cdft_tolerance_electrons",
        "minimum_snapshot_separation_ps",
        "minimum_li_contact_scale",
        "li_n_distance_angstrom",
    ):
        if not math.isfinite(float(settings[key])) or float(settings[key]) <= 0:
            raise ValueError(f"{key} must be finite and positive")
    if float(settings["cdft_tolerance_electrons"]) > 0.01:
        raise ValueError("Li cDFT residual tolerance must be <= 0.01 electron")
    for variant in settings["variants"].values():
        family = variant["xc_family"]
        fraction = float(variant.get("exact_exchange_fraction", 0.0))
        if family not in ("PBE", "PBE0") or not math.isfinite(fraction):
            raise ValueError("supported methods are PBE and PBE0")
        if (family == "PBE" and fraction != 0) or (family == "PBE0" and not 0 < fraction < 1):
            raise ValueError("exchange fraction disagrees with functional")
        for key in ("eps_scf", "cutoff_ry", "rel_cutoff_ry"):
            if not math.isfinite(float(variant[key])) or float(variant[key]) <= 0:
                raise ValueError(f"invalid method parameter: {key}")
    if settings["optimization_variant"] not in settings["variants"]:
        raise ValueError("unknown optimization method variant")
    for key in ("max_iter", "max_force", "rms_force", "max_dr", "rms_dr"):
        if not math.isfinite(float(settings["geometry"][key])) or (
            float(settings["geometry"][key]) <= 0
        ):
            raise ValueError(f"invalid geometry criterion: {key}")
    return settings


def read_cell(path: str | Path) -> np.ndarray:
    vectors = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split("#", 1)[0].split()
        if fields and fields[0].upper() in ("A", "B", "C") and len(fields) == 4:
            vectors[fields[0].upper()] = [float(x) for x in fields[1:]]
    matrix = np.asarray([vectors[key] for key in ("A", "B", "C")])
    if not np.all(np.isfinite(matrix)) or np.linalg.det(matrix) <= 0:
        raise ValueError("cell must contain finite A/B/C vectors with positive volume")
    return matrix


def seed_li(
    elements: list[str],
    positions: np.ndarray,
    cell: np.ndarray,
    donor_pairs: list[list[int]],
    *,
    distance_angstrom: float,
    minimum_scale: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Place Li near amine N donors using steric clearance, without N-H scoring."""
    if any(element not in VALENCE or element == "LI" for element in elements):
        raise ValueError("source must contain only Li-free H/C/N/O solvent")
    contact = {"H": 1.2, "C": 1.8, "N": 1.8, "O": 1.7}
    limits = np.asarray([contact[element] for element in elements])
    candidates: list[tuple[np.ndarray, list[int]]] = []
    for pair in donor_pairs:
        if len(pair) != 2 or any(elements[index] != "N" for index in pair):
            raise ValueError("bidentate donors must be pairs of nitrogen atom indices")
        first, second = pair
        delta = minimum_image_vectors(
            positions[first : first + 1], positions[second : second + 1], cell
        )[0, 0]
        length = float(np.linalg.norm(delta))
        if not 0.1 < length < 2 * distance_angstrom:
            continue
        axis = delta / length
        reference = np.eye(3)[int(np.argmin(np.abs(axis)))]
        u = np.cross(axis, reference)
        u /= np.linalg.norm(u)
        v = np.cross(axis, u)
        midpoint = positions[first] + delta / 2
        height = math.sqrt(distance_angstrom**2 - (length / 2) ** 2)
        for angle in np.linspace(0, 2 * math.pi, 36, endpoint=False):
            candidates.append(
                (midpoint + height * (math.cos(angle) * u + math.sin(angle) * v), pair)
            )

    def best(records: list[tuple[np.ndarray, list[int]]]):
        accepted = []
        for point, donors in records:
            distances = np.linalg.norm(
                minimum_image_vectors(point[None, :], positions, cell)[0], axis=1
            )
            score = float(np.min(distances / limits))
            if score >= minimum_scale:
                accepted.append((score, point, donors))
        return max(accepted, key=lambda record: record[0]) if accepted else None

    selected = best(candidates)
    mode = "bidentate"
    if selected is None:
        mode = "monodentate"
        candidates = []
        # Uniform spherical directions give both amines the same fallback rule.
        for index, element in enumerate(elements):
            if element != "N":
                continue
            for k in range(96):
                z = 1 - 2 * (k + 0.5) / 96
                angle = k * math.pi * (3 - math.sqrt(5))
                radius = math.sqrt(1 - z * z)
                direction = np.asarray([radius * math.cos(angle), radius * math.sin(angle), z])
                candidates.append((positions[index] + distance_angstrom * direction, [index]))
        selected = best(candidates)
    if selected is None:
        raise ValueError("no collision-free Li-N seed; inspect the source configuration")
    score, point, donors = selected
    point = (point @ np.linalg.inv(cell) % 1) @ cell
    return point, {
        "initial_coordination": mode,
        "donor_indices_solvent_zero_based": donors,
        "minimum_contact_scale": score,
        "selection": "maximum_steric_clearance",
    }


def render_input(
    template_path: str | Path,
    output_path: str | Path,
    *,
    coordinates: str | Path,
    cell: str | Path,
    settings: dict[str, Any],
    state: str,
    variant: str,
    optimize: bool,
) -> dict[str, Any]:
    if state not in STATES:
        raise ValueError("unknown Li electron-transfer state")
    method = settings["variants"][variant]
    elements, _positions, _ = read_xyz(coordinates)
    if elements.count("LI") != 1 or any(x not in VALENCE for x in elements):
        raise ValueError("geometry must contain one Li and H/C/N/O solvent; no ghosts")
    electrons = sum(VALENCE[x] for x in elements)
    if electrons % 2 != 1:
        raise ValueError("neutral Li/closed-shell solvent must have an odd electron count")
    matrix = read_cell(cell)
    heights = np.linalg.det(matrix) / np.asarray(
        [
            np.linalg.norm(np.cross(matrix[1], matrix[2])),
            np.linalg.norm(np.cross(matrix[2], matrix[0])),
            np.linalg.norm(np.cross(matrix[0], matrix[1])),
        ]
    )
    hybrid = method["xc_family"] == "PBE0"
    fraction = float(method.get("exact_exchange_fraction", 0))
    basis_files = "    BASIS_SET_FILE_NAME BASIS_MOLOPT"
    admm = ""
    aux = ""
    hf = ""
    if hybrid:
        if float(method["hfx_cutoff_angstrom"]) >= float(min(heights)) / 2:
            raise ValueError("HFX cutoff must be less than half the smallest cell height")
        # cFIT3 is in BASIS_ADMM, not BASIS_ADMM_MOLOPT.
        basis_files += "\n    BASIS_SET_FILE_NAME BASIS_ADMM"
        aux = f"      BASIS_SET AUX_FIT {method['aux_basis_set']}"
        admm = """    &AUXILIARY_DENSITY_MATRIX_METHOD
      METHOD BASIS_PROJECTION
      ADMM_PURIFICATION_METHOD MO_DIAG
      EXCH_CORRECTION_FUNC PBEX
    &END AUXILIARY_DENSITY_MATRIX_METHOD"""
        hf = f"""      &HF
        FRACTION {fraction:.12g}
        &SCREENING
          EPS_SCHWARZ {method["eps_schwarz"]}
          SCREEN_ON_INITIAL_P TRUE
        &END SCREENING
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE TRUNCATED
          CUTOFF_RADIUS {method["hfx_cutoff_angstrom"]}
          T_C_G_DATA t_c_g.dat
        &END INTERACTION_POTENTIAL
        &MEMORY
          MAX_MEMORY 3000
        &END MEMORY
      &END HF"""
    kinds = []
    for element in sorted(set(elements)):
        label = "Li" if element == "LI" else element
        basis = method["li_basis_set"] if element == "LI" else method["basis_set"]
        potential = "GTH-PBE-q3" if element == "LI" else f"GTH-PBE-q{VALENCE[element]}"
        kinds.append(f"""    &KIND {label}
      ELEMENT {label}
      BASIS_SET {basis}
{aux}
      POTENTIAL {potential}
    &END KIND""")
    geometry = settings["geometry"]
    motion = ""
    if optimize:
        motion = f"""&MOTION
  &GEO_OPT
    TYPE MINIMIZATION
    OPTIMIZER LBFGS
    MAX_ITER {geometry["max_iter"]}
    MAX_FORCE {geometry["max_force"]}
    RMS_FORCE {geometry["rms_force"]}
    MAX_DR {geometry["max_dr"]}
    RMS_DR {geometry["rms_dr"]}
  &END GEO_OPT
  &PRINT
    &TRAJECTORY ON
      FORMAT XYZ
      FILENAME =positions.xyz
      ADD_LAST NUMERIC
      &EACH
        GEO_OPT 1
      &END EACH
    &END TRAJECTORY
    &RESTART
      BACKUP_COPIES 2
      &EACH
        GEO_OPT 1
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION"""
    density = ""
    if not optimize:
        stride = int(settings["cube_stride"])
        density = f"""      &E_DENSITY_CUBE
        STRIDE {stride} {stride} {stride}
      &END E_DENSITY_CUBE"""
    substitutions = {
        "run_type": "GEO_OPT" if optimize else "ENERGY_FORCE",
        "basis_files": basis_files,
        "admm_block": admm,
        "hf_block": hf,
        "pbe_scale_x": f"{1 - fraction:.12g}",
        "xc_reference": method["xc_family"],
        "kinds": "\n".join(kinds),
        "motion_block": motion,
        "density_block": density,
        "coordinates_path": Path(coordinates).resolve().as_posix(),
        "cell_path": Path(cell).resolve().as_posix(),
        "li_atom_index": elements.index("LI") + 1,
        "target": settings["state_targets"][state],
        "cdft_eps_scf": settings["cdft_tolerance_electrons"],
        **method,
    }
    rendered = Template(Path(template_path).read_text()).substitute(substitutions)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered.rstrip() + "\n")
    return {
        "state": state,
        "variant": variant,
        "run_type": substitutions["run_type"],
        "charge": 0,
        "multiplicity": 2,
        "expected_valence_electrons": electrons,
        "li_atom_index": elements.index("LI") + 1,
        "target_electrons": float(settings["state_targets"][state]),
    }


def last_xyz_frame(path: str | Path) -> tuple[list[str], np.ndarray, str]:
    """Read the final complete XYZ frame; reject a truncated trailing frame."""
    lines = Path(path).read_text().splitlines()
    offset = 0
    last = None
    while offset < len(lines):
        if not lines[offset].strip():
            offset += 1
            continue
        count = int(lines[offset])
        if count < 1 or offset + count + 2 > len(lines):
            raise ValueError("truncated optimization trajectory")
        elements, positions = [], []
        for line in lines[offset + 2 : offset + count + 2]:
            fields = line.split()
            if len(fields) != 4:
                raise ValueError("invalid optimization trajectory atom record")
            elements.append(fields[0].upper())
            positions.append([float(x) for x in fields[1:]])
        last = (elements, np.asarray(positions), lines[offset + 1])
        if not np.all(np.isfinite(last[1])):
            raise ValueError("non-finite optimized coordinates")
        offset += count + 2
    if last is None:
        raise ValueError("empty optimization trajectory")
    return last


def validate_output(text: str, job: dict[str, Any], settings: dict[str, Any]) -> dict[str, Any]:
    engine = parse_cp2k_text(text)
    constraint = evaluate_cp2k_cdft_constraint(
        text,
        expected_target_electrons=float(job["target_electrons"]),
        tolerance_electrons=float(settings["cdft_tolerance_electrons"]),
    )
    problems = [*engine.problems, *constraint.problems]
    if "CDFT SCF loop FAILED" in text:
        problems.append("CP2K reported a failed CDFT outer-SCF loop")
    if engine.energy_hartree is not None and not math.isfinite(engine.energy_hartree):
        problems.append("non-finite CP2K energy")
    counts = re.findall(r"Number of electrons:\s*(\d+)", text)
    if len(counts) < 2 or [int(x) for x in counts[-2:]] != [
        (job["expected_valence_electrons"] + 1) // 2,
        (job["expected_valence_electrons"] - 1) // 2,
    ]:
        problems.append("output alpha/beta electron counts do not match the neutral doublet")
    optimized = job["run_type"] == "GEO_OPT"
    if optimized and not re.search(r"GEOMETRY OPTIMIZATION COMPLETED", text, re.I):
        problems.append("geometry optimization did not converge")
    if optimized and re.search(
        r"MAXIMUM NUMBER OF (?:OPTIMIZATION )?STEPS.*(?:REACHED|EXCEEDED)", text, re.I
    ):
        problems.append("geometry optimization exhausted its step limit")
    # CP2K's FORCE_EVAL energy includes the CDFT constraint term. Store both
    # the variational Lagrangian and the physical energy after removing it.
    # CP2K 2023.2 qs_cdft_methods.F sets energy%cdft to
    # (current_population - target_population) * constraint_strength.
    # These values are printed in each completed CDFT outer-SCF record.
    iteration = constraint.iteration
    if iteration and engine.energy_hartree is not None and (
        abs(iteration.energy_hartree - engine.energy_hartree) > 1e-7
    ):
        problems.append("final FORCE_EVAL energy does not match the final CDFT record")
    constraint_energy = (
        (iteration.current_electrons - iteration.target_electrons) * iteration.constraint_strength
        if iteration
        else None
    )
    if constraint_energy is None or not math.isfinite(constraint_energy):
        problems.append("missing/non-finite CDFT constraint energy")
    energy = engine.energy_hartree
    physical = (
        energy - constraint_energy if energy is not None and constraint_energy is not None else None
    )
    return {
        "ready": not problems,
        "problems": problems,
        "engine": engine.as_dict(),
        "constraint": constraint.as_dict(),
        "geometry_converged": optimized and not problems,
        "total_energy_hartree": energy,
        "constraint_energy_hartree": constraint_energy,
        "physical_energy_hartree": physical,
    }


def energy_metrics(energies: dict[str, float]) -> dict[str, float]:
    if set(energies) != set(PAIRS) or not all(math.isfinite(x) for x in energies.values()):
        raise ValueError("four finite matched-state energies are required")
    aa, ba, bb, ab = (energies[key] for key in PAIRS)
    convert = HARTREE_TO_EV
    return {
        "vertical_transfer_ev": (ba - aa) * convert,
        "product_relaxation_stabilization_ev": (ba - bb) * convert,
        "relaxed_transfer_ev": (bb - aa) * convert,
        "reactant_reorganization_ev": (ab - aa) * convert,
        "vertical_gap_at_product_geometry_ev": (bb - ab) * convert,
    }
