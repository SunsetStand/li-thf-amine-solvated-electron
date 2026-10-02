#!/usr/bin/env python3
"""Prepare, run and summarize the EDA/TMEDA Li relaxation comparison."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import time
import uuid
from pathlib import Path
from typing import Any

import numpy as np

from solvelec.candidates import read_xyz, write_xyz
from solvelec.composition import AVOGADRO_MOL_INV
from solvelec.li_relaxation import (
    PAIRS,
    STATUS,
    energy_metrics,
    last_xyz_frame,
    read_json,
    render_input,
    seed_li,
    settings_from,
    validate_output,
    write_json,
)
from solvelec.parsers import HARTREE_TO_EV
from solvelec.provenance import sha256_file
from solvelec.trajectory import cell_matrix, infer_element, select_representative_indices


def record(path: str | Path) -> dict[str, str]:
    file = Path(path).resolve()
    return {"path": str(file), "sha256": sha256_file(file)}


def verify(value: dict[str, str]) -> Path:
    path = Path(value["path"])
    if sha256_file(path) != value["sha256"]:
        raise ValueError(f"input checksum changed: {path}")
    return path


def source_paths(run_root: str | Path, settings: dict[str, Any], system: str) -> dict[str, Path]:
    source = settings["source_campaigns"][system]
    replica = int(settings["source_replicas"][system])
    root = Path(run_root) / source
    classical = root / "classical" / system / f"r{replica}" / "pilot"
    analysis = root / "analysis" / system / f"r{replica}"
    return {
        "analysis": analysis / "analysis.json",
        "timeseries": analysis / "timeseries.csv",
        "spec": root / "specs" / system / f"r{replica}.json",
        "tpr": classical / "production" / "production.tpr",
        "trajectory": classical / "production" / "production.xtc",
    }


def prepare(args: argparse.Namespace) -> int:
    import MDAnalysis as mda

    methods = read_json(args.methods)
    settings = settings_from(methods)
    if args.system not in settings["systems"]:
        raise ValueError("system is not selected for Li relaxation")
    paths = source_paths(args.run_root, settings, args.system)
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError(f"existing MD handoff is incomplete: {missing}; no MD rerun requested")
    analysis = read_json(paths["analysis"])
    replica = int(settings["source_replicas"][args.system])
    if (
        not analysis.get("ready")
        or analysis.get("system_id") != args.system
        or (int(analysis.get("replica", -1)) != replica)
    ):
        raise ValueError(
            "existing classical analysis must be ready for the selected system/replica"
        )
    hashes = {key: record(path) for key, path in paths.items()}
    for key in ("tpr", "trajectory", "spec"):
        if hashes[key]["sha256"] != analysis["inputs"][key]["sha256"]:
            raise ValueError(f"{key} no longer matches the accepted classical analysis")
    with paths["timeseries"].open(newline="") as handle:
        rows = [
            {key: float(value) for key, value in row.items() if value is not None}
            for row in csv.DictReader(handle)
        ]
    tau = max(
        float(value["integrated_autocorrelation_time_ps"])
        for value in analysis["metrics"]["autocorrelation"].values()
    )
    separation = max(
        float(settings["minimum_snapshot_separation_ps"]),
        float(methods["trajectory_analysis"]["decorrelation_multiplier"]) * tau,
    )
    selection = select_representative_indices(
        [{**row, "time_ps": row["elapsed_ps"]} for row in rows],
        ["density_g_ml", "volume_nm3"],
        count=int(settings["snapshots_per_system"]),
        minimum_time_ps=float(methods["trajectory_analysis"]["equilibrated_start_ns"]) * 1000,
        minimum_separation_ps=separation,
    )
    universe = mda.Universe(str(paths["tpr"]), str(paths["trajectory"]))
    elements = []
    for atom in universe.atoms:
        try:
            explicit = atom.element
        except (AttributeError, ValueError):
            explicit = None
        elements.append(infer_element(str(atom.name), explicit))
    donor_pairs = []
    for residue in universe.residues:
        indices = [int(atom.index) for atom in residue.atoms if elements[int(atom.index)] == "N"]
        if indices:
            if len(indices) != 2:
                raise ValueError("EDA/TMEDA source residues must contain two N atoms")
            donor_pairs.append(indices)
    if not donor_pairs:
        raise ValueError("no amine nitrogen pairs in the source topology")
    spec = read_json(paths["spec"])
    amine_count = int(spec["component_counts"][str(spec["amine"])])
    if len(donor_pairs) != amine_count:
        raise ValueError("source amine molecule count disagrees with the topology")
    output = Path(args.output_dir).resolve()
    seeds = []
    for ordinal, selected in enumerate(selection, 1):
        row = rows[selected]
        frame_index = int(row["frame_index"])
        frame = universe.trajectory[frame_index]
        if abs(float(frame.time) - row["time_ps"]) > 0.01:
            raise ValueError("selected trajectory frame time disagrees with analysis")
        matrix = cell_matrix(frame.dimensions)
        # Retain whole molecules when possible; quantum calculations are periodic.
        try:
            positions = np.asarray(universe.atoms.unwrap(compound="residues", inplace=False))
        except (AttributeError, ValueError):
            positions = np.asarray(universe.atoms.positions, dtype=float).copy()
        li, placement = seed_li(
            elements,
            positions,
            matrix,
            donor_pairs,
            distance_angstrom=float(settings["li_n_distance_angstrom"]),
            minimum_scale=float(settings["minimum_li_contact_scale"]),
        )
        snapshot = f"s{ordinal:02d}"
        directory = output / snapshot
        xyz = directory / "coordinates.xyz"
        cell = directory / "cell.inc"
        write_xyz(
            xyz,
            ["LI", *elements],
            np.vstack([li, positions]),
            f"system={args.system} snapshot={snapshot} frame={frame_index}",
        )
        vectors = ["&CELL"]
        for label, vector in zip(("A", "B", "C"), matrix, strict=True):
            vectors.append(f"  {label} " + " ".join(f"{x:.10f}" for x in vector))
        vectors.extend(["  PERIODIC XYZ", "&END CELL"])
        cell.write_text("\n".join(vectors) + "\n")
        volume_nm3 = float(np.linalg.det(matrix)) / 1000
        seeds.append(
            {
                "snapshot": snapshot,
                "frame_index": frame_index,
                "time_ps": row["time_ps"],
                "elapsed_ps": row["elapsed_ps"],
                "structure": {
                    "xyz": record(xyz),
                    "cell": record(cell),
                    "cell_vectors_angstrom": matrix.tolist(),
                },
                "placement": placement,
                "amine_count": amine_count,
                "volume_nm3": volume_nm3,
                "achieved_amine_concentration_m": amine_count
                / (AVOGADRO_MOL_INV * volume_nm3 * 1e-24),
                "li_concentration_m": 1 / (AVOGADRO_MOL_INV * volume_nm3 * 1e-24),
            }
        )
    write_json(
        output / "bank.json",
        {
            "schema_version": 1,
            "ready": True,
            "scientific_status": STATUS,
            "system_id": args.system,
            "source_campaign": settings["source_campaigns"][args.system],
            "source_replica": replica,
            "snapshot_count": len(seeds),
            "seeds": seeds,
            "minimum_separation_ps": separation,
            "source": hashes,
            "methods": record(args.methods),
            "component_counts": spec["component_counts"],
        },
    )
    return 0


def render(args: argparse.Namespace) -> int:
    settings = settings_from(read_json(args.methods))
    bank = read_json(args.bank)
    matches = [seed for seed in bank["seeds"] if seed["snapshot"] == args.snapshot]
    if not bank.get("ready") or len(matches) != 1:
        raise ValueError("selected snapshot is missing or duplicated")
    seed = matches[0]
    verify(seed["structure"]["cell"])
    reference_elements, _, _ = read_xyz(verify(seed["structure"]["xyz"]))
    if args.optimize:
        expected_geometry = "seed" if args.state == "a" else "a"
        if args.geometry_state != expected_geometry:
            raise ValueError("optimization must follow seed -> state a -> state b")
    elif args.geometry_state not in ("a", "b"):
        raise ValueError("single points require an optimized a or b geometry")
    if args.geometry_state == "seed":
        if sha256_file(args.coordinates) != seed["structure"]["xyz"]["sha256"]:
            raise ValueError("initial optimization coordinates do not match the selected seed")
    elif not args.geometry_validation:
        raise ValueError("optimized geometry validation is required")
    elements, _, _ = read_xyz(args.coordinates)
    if elements != reference_elements:
        raise ValueError("optimization changed the atom identities or ordering")
    if args.geometry_validation:
        validation = read_json(args.geometry_validation)
        if not validation.get("ready"):
            raise ValueError("upstream geometry optimization is not ready")
        if validation["job"]["system_id"] != bank["system_id"] or (
            validation["job"]["snapshot"] != args.snapshot
        ):
            raise ValueError("geometry belongs to a different system/snapshot")
        if validation["job"]["state"] != args.geometry_state:
            raise ValueError("geometry state does not match its optimization")
        if validation["job"]["variant"] != settings["optimization_variant"]:
            raise ValueError("geometry was optimized with a different method")
        if validation["job"]["methods"]["sha256"] != sha256_file(args.methods):
            raise ValueError("geometry was optimized with a different configuration")
        if validation["optimized_xyz"]["sha256"] != sha256_file(args.coordinates):
            raise ValueError("optimized coordinates do not match their validation")
        if validation["job"]["cell"]["sha256"] != seed["structure"]["cell"]["sha256"]:
            raise ValueError("upstream optimization changed the cell")
    job = render_input(
        args.template,
        args.output,
        coordinates=args.coordinates,
        cell=seed["structure"]["cell"]["path"],
        settings=settings,
        state=args.state,
        variant=args.variant,
        optimize=args.optimize,
    )
    job.update(
        {
            "schema_version": 1,
            "system_id": bank["system_id"],
            "snapshot": args.snapshot,
            "geometry_state": args.geometry_state,
            "coordinates": record(args.coordinates),
            "cell": seed["structure"]["cell"],
            "methods": record(args.methods),
            "input": record(args.output),
            "bank": record(args.bank),
        }
    )
    if args.geometry_validation:
        job["geometry_validation"] = record(args.geometry_validation)
    write_json(Path(args.output).with_suffix(".job.json"), job)
    return 0


def check_data(input_text: str) -> None:
    import re

    data_dir = Path(os.environ["CP2K_DATA_DIR"])
    for match in re.finditer(
        r"(?:BASIS_SET_FILE_NAME|POTENTIAL_FILE_NAME|PARAMETER_FILE_NAME|T_C_G_DATA)\s+(\S+)",
        input_text,
    ):
        if not (data_dir / match.group(1)).is_file():
            raise ValueError(f"CP2K data file is unavailable: {data_dir / match.group(1)}")
    orbital = (data_dir / "BASIS_MOLOPT").read_text()
    auxiliary = (data_dir / "BASIS_ADMM").read_text() if "&HF" in input_text else ""
    for match in re.finditer(r"&KIND\s+(\S+)(.*?)&END KIND", input_text, re.S):
        element, body = match.groups()
        for line in body.splitlines():
            fields = line.split()
            if fields[:1] != ["BASIS_SET"]:
                continue
            basis = fields[-1]
            data = auxiliary if "AUX_FIT" in fields else orbital
            if not any(
                fields[0].upper() == element.upper() and basis in fields[1:]
                for row in data.splitlines()
                if (fields := row.split())
            ):
                raise ValueError(
                    f"CP2K basis is absent from the selected data file: {element} {basis}"
                )


def run(args: argparse.Namespace) -> int:
    if not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("CP2K execution requires a Slurm allocation")
    job = read_json(args.job)
    for key in ("coordinates", "cell", "methods", "input", "bank"):
        verify(job[key])
    if "geometry_validation" in job:
        verify(job["geometry_validation"])
    settings = settings_from(read_json(job["methods"]["path"]))
    command = args.engine_command[1:] if args.engine_command[:1] == ["--"] else args.engine_command
    if not command:
        raise ValueError("engine launcher required after --")
    directory = Path(args.workdir).resolve()
    attempt = directory / "attempts" / f"{time.time_ns()}-{uuid.uuid4().hex[:8]}"
    attempt.mkdir(parents=True)
    input_path = attempt / "cp2k.inp"
    shutil.copy2(job["input"]["path"], input_path)
    check_data(input_path.read_text())
    raw = attempt / "cp2k.out"
    environment = {**os.environ, "OMP_NUM_THREADS": "1"}
    with raw.open("w") as handle:
        completed = subprocess.run(
            [*command, "-i", str(input_path)],
            cwd=attempt,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    result = validate_output(raw.read_text(errors="replace"), job, settings)
    if completed.returncode:
        result["problems"].append(f"engine exit code {completed.returncode}")
    result.update({"job": job, "raw_output": record(raw), "attempt_directory": str(attempt)})
    optimized = job["run_type"] == "GEO_OPT"
    if not result["problems"] and optimized:
        try:
            elements, positions, comment = last_xyz_frame(attempt / "positions.xyz")
            expected, _, _ = read_xyz(job["coordinates"]["path"])
            if elements != expected:
                raise ValueError("final trajectory atom identities/order changed")
            # XYZ comments print total FORCE_EVAL energy in Hartree.
            import re

            match = re.search(r"\bE\s*=\s*([-+0-9.Ee]+)", comment)
            if match is None or abs(float(match.group(1)) - result["total_energy_hartree"]) > 1e-7:
                raise ValueError("final trajectory frame does not match the final evaluated energy")
            write_xyz(directory / "optimized.xyz", elements, positions, comment)
            result["optimized_xyz"] = record(directory / "optimized.xyz")
        except (OSError, ValueError) as exc:
            result["problems"].append(str(exc))
    if not result["problems"] and not optimized:
        for label, marker in (("electron", "ELECTRON_DENSITY"), ("spin", "SPIN_DENSITY")):
            cubes = sorted(attempt.glob(f"*{marker}*.cube"))
            if len(cubes) != 1:
                result["problems"].append(f"expected one {label} cube, found {len(cubes)}")
                continue
            shutil.copy2(cubes[0], directory / f"{label}.cube")
            result[f"{label}_cube"] = record(directory / f"{label}.cube")
    result["ready"] = not result["problems"]
    write_json(attempt / "validation.json", result)
    if not result["ready"]:
        shutil.copy2(raw, directory / "cp2k.out.failed")
        write_json(directory / "validation.failed.json", result)
        print(f"CP2K validation failed: {result['problems']}; logs: {attempt}", file=sys.stderr)
        return completed.returncode or 4
    shutil.copy2(raw, directory / "cp2k.out")
    write_json(directory / "validation.json", result)
    return 0


def analyze(args: argparse.Namespace) -> int:
    settings = settings_from(read_json(args.methods))
    bank = read_json(args.bank)
    root = Path(args.snapshot_dir)
    problems = []
    energies = {}
    provenance = {}
    for pair in PAIRS:
        path = root / "sp" / args.variant / pair / "validation.json"
        value = read_json(path)
        job = value["job"]
        if not value.get("ready"):
            problems.append(f"{pair}: CP2K calculation is not ready")
        if (
            job["system_id"],
            job["snapshot"],
            job["variant"],
            job["state"],
            job["geometry_state"],
        ) != (bank["system_id"], args.snapshot, args.variant, pair[0], pair[1]):
            raise ValueError(f"{pair}: mixed electronic states, geometries or methods")
        verify(job["input"])
        verify(job["coordinates"])
        verify(job["cell"])
        verify(job["bank"])
        verify(job["geometry_validation"])
        if job["bank"]["sha256"] != sha256_file(args.bank):
            raise ValueError("single points used a different snapshot bank")
        if job["methods"]["sha256"] != sha256_file(args.methods):
            raise ValueError("single points used a different methods configuration")
        verify(value["raw_output"])
        checked = validate_output(Path(value["raw_output"]["path"]).read_text(), job, settings)
        if not checked["ready"]:
            problems.extend(checked["problems"])
        energies[pair] = checked["physical_energy_hartree"]
        provenance[pair] = {"validation": record(path), "job": job, "output": value["raw_output"]}
    for first, second in (("aa", "ba"), ("bb", "ab")):
        if (
            provenance[first]["job"]["coordinates"]["sha256"]
            != (provenance[second]["job"]["coordinates"]["sha256"])
        ):
            raise ValueError("cross-state single points are not on identical coordinates")
    if len({provenance[key]["job"]["cell"]["sha256"] for key in PAIRS}) != 1:
        raise ValueError("energy matrix mixes cells")
    metrics = energy_metrics(energies) if not problems else None
    diagnostics = []
    if metrics and args.variant == settings["optimization_variant"]:
        for state in ("a", "b"):
            pair = state + state
            validation = read_json(verify(provenance[pair]["job"]["geometry_validation"]))
            optimized_energy = validation["physical_energy_hartree"]
            difference = abs(energies[pair] - optimized_energy) * HARTREE_TO_EV
            if difference > float(settings["energy_tolerance_ev"]):
                problems.append(f"{pair}: single-point energy differs from its optimized state")
        if metrics["product_relaxation_stabilization_ev"] < -float(settings["energy_tolerance_ev"]):
            problems.append("product relaxation increased the energy from its starting geometry")
        if metrics["reactant_reorganization_ev"] < -float(settings["energy_tolerance_ev"]):
            diagnostics.append(
                "product geometry has a lower reference-state energy than its local minimum"
            )
    write_json(
        args.output,
        {
            "schema_version": 1,
            "ready": not problems,
            "scientific_status": STATUS,
            "system_id": bank["system_id"],
            "snapshot": args.snapshot,
            "variant": args.variant,
            "geometry_variant": settings["optimization_variant"],
            "physical_energies_hartree": energies,
            "metrics_ev": metrics,
            "problems": problems,
            "diagnostics": diagnostics,
            "provenance": provenance,
            "bank": record(args.bank),
        },
    )
    return 0


def summary(args: argparse.Namespace) -> int:
    settings = settings_from(read_json(args.methods))
    expected = {
        (system, snapshot, variant)
        for system in settings["systems"]
        for snapshot in args.snapshots
        for variant in args.variants
    }
    records = [read_json(path) for path in args.records]
    keys = [(value["system_id"], value["snapshot"], value["variant"]) for value in records]
    if len(keys) != len(set(keys)) or set(keys) != expected:
        raise ValueError("summary must contain each requested system/snapshot/method exactly once")
    means = {}
    for variant in args.variants:
        means[variant] = {}
        for system in settings["systems"]:
            group = [
                value
                for value in records
                if value["system_id"] == system and value["variant"] == variant
            ]
            ready = all(value.get("ready") for value in group)
            descriptors = {}
            if ready:
                for name in group[0]["metrics_ev"]:
                    values = [float(value["metrics_ev"][name]) for value in group]
                    descriptors[name] = {
                        "values_ev": values,
                        "mean_ev": float(np.mean(values)),
                        "standard_deviation_ev": float(np.std(values, ddof=1))
                        if len(values) > 1
                        else None,
                    }
            means[variant][system] = {
                "ready": ready,
                "snapshot_count": len(group),
                "metrics": descriptors,
            }
    comparisons = {}
    for variant, systems in means.items():
        if all(value["ready"] for value in systems.values()):
            comparisons[variant] = {
                name: systems["eda_3m"]["metrics"][name]["mean_ev"]
                - systems["tmeda_3m"]["metrics"][name]["mean_ev"]
                for name in systems["eda_3m"]["metrics"]
            }
    ready = all(value.get("ready") for value in records)
    write_json(
        args.output,
        {
            "schema_version": 1,
            "ready": ready,
            "campaign": args.campaign,
            "scientific_status": STATUS,
            "kind": "li_relaxation",
            "geometry_variant": settings["optimization_variant"],
            "record_count": len(records),
            "systems_by_variant": means,
            "eda_minus_tmeda_ev": comparisons,
            "records": records,
            "energy_definition": (
                "Physical DFT energy, with residual CDFT constraint energy removed"
            ),
            "interpretation": (
                "Fixed-cell local-minimum reaction energies; sample spread is reported."
            ),
            "methods": record(args.methods),
        },
    )
    csv_path = Path(args.output).with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    names = sorted({name for value in records for name in (value["metrics_ev"] or {})})
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["system", "snapshot", "variant", "ready", *names]
        )
        writer.writeheader()
        for value in records:
            writer.writerow(
                {
                    "system": value["system_id"],
                    "snapshot": value["snapshot"],
                    "variant": value["variant"],
                    "ready": value["ready"],
                    **(value["metrics_ev"] or {}),
                }
            )
    return 0


def gate(args: argparse.Namespace) -> int:
    if not read_json(args.summary).get("ready"):
        raise ValueError("Li relaxation calculations did not all pass")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        f"sha256 {hashlib.sha256(Path(args.summary).read_bytes()).hexdigest()}  "
        f"{Path(args.summary).name}\n"
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    p = commands.add_parser("prepare")
    for name in ("run-root", "methods", "system", "output-dir"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(func=prepare)
    p = commands.add_parser("render")
    for name in (
        "methods",
        "bank",
        "snapshot",
        "state",
        "variant",
        "coordinates",
        "template",
        "output",
    ):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--geometry-state", default="seed")
    p.add_argument("--geometry-validation")
    p.add_argument("--optimize", action="store_true")
    p.set_defaults(func=render)
    p = commands.add_parser("run")
    p.add_argument("--job", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("engine_command", nargs=argparse.REMAINDER)
    p.set_defaults(func=run)
    p = commands.add_parser("analyze")
    for name in ("methods", "bank", "snapshot", "variant", "snapshot-dir", "output"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(func=analyze)
    p = commands.add_parser("summary")
    for name in ("campaign", "methods", "output"):
        p.add_argument(f"--{name}", required=True)
    for name in ("snapshots", "variants", "records"):
        p.add_argument(f"--{name}", required=True, nargs="+")
    p.set_defaults(func=summary)
    p = commands.add_parser("gate")
    p.add_argument("--summary", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=gate)
    args = parser.parse_args()
    try:
        return int(args.func(args))
    except (KeyError, OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
