from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, cast

from .classification import LocalizationMetrics, classify_localization
from .composition import AVOGADRO_MOL_INV, NM3_TO_L, concentration_molar
from .config import (
    campaign_matrix,
    load_repository_configs,
    make_system_spec,
    repository_root,
    validate_repository_configs,
)
from .cube import analyze_spin_density, read_cube
from .engines import doctor_report
from .parsers import parse_output
from .provenance import write_manifest
from .rendering import render_cp2k, render_orca, render_packmol


def _json_dump(value: Any) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def _root(value: str | None) -> Path:
    return Path(value).resolve() if value else repository_root()


def cmd_validate(args: argparse.Namespace) -> int:
    errors = validate_repository_configs(_root(args.root))
    result = {"valid": not errors, "errors": errors}
    _json_dump(result)
    return 0 if not errors else 2


def cmd_doctor(args: argparse.Namespace) -> int:
    report = doctor_report()
    _json_dump(report)
    if args.strict_engines:
        engines = cast(list[dict[str, object]], report["engines"])
        chemistry_missing = [
            item["name"]
            for item in engines
            if item["category"] == "chemistry" and not item["found"]
        ]
        return 3 if chemistry_missing else 0
    return 0


def _matrix_records(root: Path, campaign_name: str) -> list[dict[str, Any]]:
    campaign, systems, _ = load_repository_configs(root)
    records = []
    for spec, replica in campaign_matrix(campaign_name, campaign, systems):
        record = spec.as_dict()
        record["replica"] = replica
        records.append(record)
    return records


def cmd_matrix(args: argparse.Namespace) -> int:
    _json_dump(
        {"campaign": args.campaign, "jobs": _matrix_records(_root(args.root), args.campaign)}
    )
    return 0


def cmd_write_spec(args: argparse.Namespace) -> int:
    root = _root(args.root)
    campaign, systems, _ = load_repository_configs(root)
    spec = make_system_spec(args.system, campaign, systems).as_dict()
    spec["replica"] = args.replica
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(spec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


def _initial_box_angstrom(spec: dict[str, Any], systems: dict[str, Any]) -> float:
    volume_l_mol_scaled = spec["thf_count"] * float(systems["thf"]["molar_volume_l_mol"])
    if spec["amine"]:
        volume_l_mol_scaled += spec["amine_count_initial"] * float(
            systems["amines"][spec["amine"]]["molar_volume_l_mol"]
        )
    volume_l = volume_l_mol_scaled / AVOGADRO_MOL_INV
    volume_nm3 = volume_l / NM3_TO_L
    return volume_nm3 ** (1.0 / 3.0) * 10.0


def cmd_render_packmol(args: argparse.Namespace) -> int:
    root = _root(args.root)
    campaign, systems, _ = load_repository_configs(root)
    spec = make_system_spec(args.system, campaign, systems).as_dict()
    box = args.box_angstrom or _initial_box_angstrom(spec, systems)
    amine_structure = None
    if spec["amine"]:
        amine_structure = f"molecules/{spec['amine']}.pdb"
    render_packmol(
        args.output,
        output_pdb=f"{args.system}_r{args.replica}.pdb",
        box_angstrom=box,
        thf_structure="molecules/thf.pdb",
        thf_count=int(spec["thf_count"]),
        amine_structure=amine_structure,
        amine_count=int(spec["amine_count_initial"]),
        seed=args.seed or args.replica,
    )
    return 0


def cmd_render_cp2k(args: argparse.Namespace) -> int:
    root = _root(args.root)
    _, _, methods = load_repository_configs(root)
    render_cp2k(
        root / "workflow" / "templates" / "cp2k" / "pbe0_cdft.inp.tpl",
        args.output,
        args.state,
        args.project,
        args.coordinates_include,
        args.cell_include,
        methods["cp2k"],
        args.li_atom_index,
        args.constrained,
    )
    return 0


def cmd_render_orca(args: argparse.Namespace) -> int:
    root = _root(args.root)
    _, _, methods = load_repository_configs(root)
    coordinates = Path(args.coordinates).read_text(encoding="utf-8")
    lines = coordinates.splitlines()
    if lines and lines[0].strip().isdigit():
        coordinates = "\n".join(lines[2:])
    render_orca(
        root / "workflow" / "templates" / "orca" / "delta_scf.inp.tpl",
        args.output,
        args.state,
        coordinates,
        methods["orca"],
    )
    return 0


def cmd_analyze_cube(args: argparse.Namespace) -> int:
    metrics = analyze_spin_density(read_cube(args.cube), clip_negative=not args.absolute)
    _json_dump(metrics.as_dict())
    return 0


def cmd_classify(args: argparse.Namespace) -> int:
    root = _root(args.root)
    _, _, methods = load_repository_configs(root)
    metrics = LocalizationMetrics(
        electron_count=args.electron_count,
        li_spin=args.li_spin,
        max_atomic_spin=args.max_atomic_spin,
        interstitial_fraction=args.interstitial_fraction,
    )
    _json_dump(classify_localization(metrics, methods["localization_thresholds"]).as_dict())
    return 0


def cmd_parse_output(args: argparse.Namespace) -> int:
    result = parse_output(args.engine, args.path)
    _json_dump(result.as_dict())
    return 0 if result.converged else 4


def cmd_manifest(args: argparse.Namespace) -> int:
    write_manifest(args.output, _root(args.root), args.inputs, args.campaign)
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    root = _root(args.root)
    records = _matrix_records(root, args.campaign)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    systems = sorted({str(record["system_id"]) for record in records})
    text = [
        f"# Campaign report: {args.campaign}",
        "",
        "> This is an input/provenance readiness report. "
        "It contains no simulated scientific result.",
        "",
        f"- Systems: {len(systems)}",
        f"- Planned replica jobs: {len(records)}",
        "- Status: inputs planned; chemistry-engine gates remain external",
        "",
        "## Systems",
        "",
        *[f"- `{system}`" for system in systems],
        "",
    ]
    output.write_text("\n".join(text), encoding="utf-8")
    return 0


def cmd_test(args: argparse.Namespace) -> int:
    root = _root(args.root)
    command = [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"]
    return subprocess.run(command, cwd=root, check=False).returncode


def cmd_concentration(args: argparse.Namespace) -> int:
    _json_dump(
        {
            "molecule_count": args.count,
            "volume_nm3": args.volume_nm3,
            "concentration_m": concentration_molar(args.count, args.volume_nm3),
        }
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="solvelec")
    parser.add_argument("--root", help="repository root (auto-detected by default)")
    sub = parser.add_subparsers(dest="command", required=True)

    validate = sub.add_parser("validate", help="validate committed configuration and safety gates")
    validate.set_defaults(func=cmd_validate)

    doctor = sub.add_parser("doctor", help="report workflow and chemistry-engine capabilities")
    doctor.add_argument("--strict-engines", action="store_true")
    doctor.set_defaults(func=cmd_doctor)

    matrix = sub.add_parser("matrix", help="emit the expanded campaign matrix")
    matrix.add_argument("--campaign", default="pilot")
    matrix.set_defaults(func=cmd_matrix)

    write_spec = sub.add_parser("write-spec")
    write_spec.add_argument("--system", required=True)
    write_spec.add_argument("--replica", type=int, required=True)
    write_spec.add_argument("--output", required=True)
    write_spec.set_defaults(func=cmd_write_spec)

    packmol = sub.add_parser("render-packmol")
    packmol.add_argument("--system", required=True)
    packmol.add_argument("--replica", type=int, required=True)
    packmol.add_argument("--output", required=True)
    packmol.add_argument("--box-angstrom", type=float)
    packmol.add_argument("--seed", type=int)
    packmol.set_defaults(func=cmd_render_packmol)

    cp2k = sub.add_parser("render-cp2k")
    cp2k.add_argument("--state", choices=["solvated_electron", "detached"], required=True)
    cp2k.add_argument("--project", required=True)
    cp2k.add_argument("--coordinates-include", required=True)
    cp2k.add_argument("--cell-include", required=True)
    cp2k.add_argument("--li-atom-index", type=int, required=True)
    cp2k.add_argument("--constrained", action="store_true")
    cp2k.add_argument("--output", required=True)
    cp2k.set_defaults(func=cmd_render_cp2k)

    orca = sub.add_parser("render-orca")
    orca.add_argument("--state", choices=["solvated_electron", "detached"], required=True)
    orca.add_argument("--coordinates", required=True)
    orca.add_argument("--output", required=True)
    orca.set_defaults(func=cmd_render_orca)

    cube = sub.add_parser("analyze-cube")
    cube.add_argument("cube")
    cube.add_argument(
        "--absolute", action="store_true", help="use absolute rather than positive density"
    )
    cube.set_defaults(func=cmd_analyze_cube)

    classify = sub.add_parser("classify")
    classify.add_argument("--electron-count", type=float, required=True)
    classify.add_argument("--li-spin", type=float, required=True)
    classify.add_argument("--max-atomic-spin", type=float, required=True)
    classify.add_argument("--interstitial-fraction", type=float, required=True)
    classify.set_defaults(func=cmd_classify)

    output_parser = sub.add_parser("parse-output")
    output_parser.add_argument("--engine", choices=["cp2k", "orca"], required=True)
    output_parser.add_argument("path")
    output_parser.set_defaults(func=cmd_parse_output)

    manifest = sub.add_parser("manifest")
    manifest.add_argument("--output", required=True)
    manifest.add_argument("--campaign")
    manifest.add_argument("inputs", nargs="*")
    manifest.set_defaults(func=cmd_manifest)

    report = sub.add_parser("report")
    report.add_argument("--campaign", default="pilot")
    report.add_argument("--output", required=True)
    report.set_defaults(func=cmd_report)

    tests = sub.add_parser("test")
    tests.set_defaults(func=cmd_test)

    concentration = sub.add_parser("concentration")
    concentration.add_argument("--count", type=int, required=True)
    concentration.add_argument("--volume-nm3", type=float, required=True)
    concentration.set_defaults(func=cmd_concentration)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (ValueError, KeyError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
