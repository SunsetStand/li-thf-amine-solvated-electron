#!/usr/bin/env python3
"""Benchmark Stage-B2C rich/poor ordering against tighter PBE and PBE0."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

from solvelec.candidates import read_xyz
from solvelec.provenance import sha256_file
from solvelec.rendering import render_stage_b2c_method_benchmark_cp2k

SCIENTIFIC_STATUS = "NUMERICAL_METHOD_SENSITIVITY_BENCHMARK_ONLY"


def _load_base_module():
    path = Path(__file__).with_name("prepare_stage_b2c_preferential.py")
    spec = importlib.util.spec_from_file_location("stage_b2c_preferential_base", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_json(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def _write_json(path: str | Path, value: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)


def _validate_gate(summary_path: str | Path, gate_path: str | Path) -> dict[str, Any]:
    summary = Path(summary_path)
    gate = Path(gate_path)
    if not summary.is_file() or not gate.is_file():
        raise ValueError("method benchmark requires an accepted production summary and gate")
    value = _read_json(summary)
    if not value.get("ready"):
        raise ValueError("method-benchmark source summary is not ready")
    digest = hashlib.sha256(summary.read_bytes()).hexdigest()
    if gate.read_text(encoding="utf-8").split() != ["sha256", digest, summary.name]:
        raise ValueError("method-benchmark source gate does not match its summary")
    return value


def _settings(methods: dict[str, Any]) -> dict[str, Any]:
    settings = methods["stage_b2c_method_benchmark"]
    if settings.get("scientific_status") != SCIENTIFIC_STATUS:
        raise ValueError("method benchmark is not explicitly diagnostic-only")
    return settings


def _variant(methods: dict[str, Any], variant_id: str) -> dict[str, Any]:
    benchmark = _settings(methods)
    records = [record for record in benchmark["variants"] if record.get("id") == variant_id]
    if len(records) != 1:
        raise ValueError(f"unknown or duplicate method variant: {variant_id}")
    return {
        **methods["stage_b2c_preferential_smoke"],
        **records[0],
        "scientific_status": SCIENTIFIC_STATUS,
    }


def _preference(delta_ev: float | None, tolerance_ev: float) -> str:
    if delta_ev is None:
        return "unavailable"
    if abs(delta_ev) <= tolerance_ev:
        return "indistinguishable_at_benchmark_tolerance"
    return "eda_rich" if delta_ev > 0 else "eda_poor"


def run_plan(args: argparse.Namespace) -> int:
    methods_path = Path(args.methods).resolve()
    methods = _read_json(methods_path)
    settings = _settings(methods)
    production = _validate_gate(args.production_summary, args.production_gate)
    if production.get("scientific_status") != settings["source_scientific_status"]:
        raise ValueError("method-benchmark source has the wrong scientific status")
    comparisons = list(production.get("comparisons", []))
    selected: list[dict[str, Any]] = []
    for system in settings["systems"]:
        group = [record for record in comparisons if record.get("system_id") == system]
        if not group or not all(record.get("ready") for record in group):
            raise ValueError(f"production comparison bank is incomplete for {system}")
        if any(
            record.get("vertical_attachment_proxy_ev", {}).get("rich_minus_poor") is None
            for record in group
        ):
            raise ValueError(f"production comparison bank lacks an energy delta for {system}")
        maximum = max(
            group,
            key=lambda record: (
                float(record["vertical_attachment_proxy_ev"]["rich_minus_poor"]),
                -int(record["replica"]),
            ),
        )
        minimum = min(
            group,
            key=lambda record: (
                float(record["vertical_attachment_proxy_ev"]["rich_minus_poor"]),
                int(record["replica"]),
            ),
        )
        extrema = [maximum, minimum]
        replicas = [int(record["replica"]) for record in extrema]
        configured = [int(value) for value in settings["representative_replicas"][system]]
        if len(set(replicas)) != 2:
            raise ValueError(f"positive/negative representatives are not distinct for {system}")
        if replicas != configured:
            raise ValueError(
                f"configured representatives for {system} are stale: {configured} != {replicas}"
            )
        for label, record in zip(settings["selection_labels"], extrema, strict=True):
            selected.append(
                {
                    "system_id": system,
                    "replica": int(record["replica"]),
                    "selection_label": str(label),
                    "baseline_comparison": record,
                }
            )
    variants = [str(record["id"]) for record in settings["variants"]]
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "ready": True,
            "campaign": args.campaign,
            "kind": "stage_b2c_method_benchmark_plan",
            "scientific_status": SCIENTIFIC_STATUS,
            "selected_pairs": selected,
            "method_variants": variants,
            "expected_pair_analyses": len(selected) * len(variants) * 2,
            "expected_cp2k_single_points": len(selected) * len(variants) * 2 * 2,
            "source": {
                "production_summary": {
                    "path": str(Path(args.production_summary).resolve()),
                    "sha256": sha256_file(args.production_summary),
                },
                "production_gate": str(Path(args.production_gate).resolve()),
                "methods": {"path": str(methods_path), "sha256": sha256_file(methods_path)},
            },
            "interpretation": settings["interpretation"],
        },
    )
    return 0


def _selected_pair(plan: dict[str, Any], system: str, replica: int) -> dict[str, Any]:
    matches = [
        record
        for record in plan.get("selected_pairs", [])
        if record.get("system_id") == system and int(record.get("replica", -1)) == replica
    ]
    if not plan.get("ready") or len(matches) != 1:
        raise ValueError(f"{system} r{replica} is not uniquely selected by the benchmark plan")
    return matches[0]


def run_render(args: argparse.Namespace) -> int:
    base = _load_base_module()
    methods = _read_json(args.methods)
    settings = _variant(methods, args.method_variant)
    plan = _read_json(args.plan)
    _selected_pair(plan, args.system, int(args.replica))
    manifest = _read_json(args.manifest)
    seed = base._ready_seed(manifest, args.seed_role)
    structure = seed["structure"]
    coordinates = Path(structure["xyz"]["path"])
    cell = Path(structure["cell"]["path"])
    if sha256_file(coordinates) != structure["xyz"]["sha256"]:
        raise ValueError("method-benchmark coordinate checksum mismatch")
    if sha256_file(cell) != structure["cell"]["sha256"]:
        raise ValueError("method-benchmark cell checksum mismatch")
    elements, _positions, _comment = read_xyz(coordinates)
    if any(element.upper() == "LI" for element in elements):
        raise ValueError("method-benchmark coordinates unexpectedly contain Li")
    states = base._states(methods["stage_b2c_preferential_smoke"])
    for state_id, output in (("neutral", args.neutral_output), ("anion", args.anion_output)):
        render_stage_b2c_method_benchmark_cp2k(
            args.template,
            output,
            project=f"{args.project_prefix}_{state_id}_stage_b2c_method_benchmark",
            coordinates_path=coordinates,
            cell_path=cell,
            method=settings,
            state=states[state_id],
        )
    return 0


def _method_signature_problems(text: str, method: dict[str, Any]) -> list[str]:
    upper = text.upper()
    problems: list[str] = []
    if str(method["basis_set"]).upper() not in upper:
        problems.append("configured orbital basis is absent from CP2K input")
    if f"CUTOFF {method['cutoff_ry']}" not in upper:
        problems.append("configured plane-wave cutoff is absent from CP2K input")
    if method["xc_family"] == "PBE0":
        for marker in ("&HF", "AUXILIARY_DENSITY_MATRIX_METHOD", "BASIS_ADMM_MOLOPT"):
            if marker not in upper:
                problems.append(f"PBE0 input lacks {marker}")
    elif "&HF" in upper or "AUXILIARY_DENSITY_MATRIX_METHOD" in upper:
        problems.append("PBE input unexpectedly enables hybrid/ADMM machinery")
    return problems


def run_analyze(args: argparse.Namespace) -> int:
    base = _load_base_module()
    methods = _read_json(args.methods)
    method = _variant(methods, args.method_variant)
    plan = _read_json(args.plan)
    selection = _selected_pair(plan, args.system, int(args.replica))
    result = int(base.run_analyze(args))
    if result:
        return result
    output = _read_json(args.output)
    problems = list(output.get("problems", []))
    for path in (args.neutral_input, args.anion_input):
        problems.extend(
            _method_signature_problems(Path(path).read_text(encoding="utf-8"), method)
        )
    output.update(
        {
            "kind": "stage_b2c_method_benchmark_pair",
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "method_variant": args.method_variant,
            "selection_label": selection["selection_label"],
            "problems": problems,
            "method": {
                key: method.get(key)
                for key in (
                    "id",
                    "functional",
                    "xc_family",
                    "basis_set",
                    "ghost_basis_set",
                    "aux_basis_set",
                    "cutoff_ry",
                    "rel_cutoff_ry",
                    "eps_scf",
                    "max_scf",
                    "admm",
                    "exact_exchange_fraction",
                    "hfx_cutoff_angstrom",
                )
                if key in method
            },
            "benchmark_plan": {
                "path": str(Path(args.plan).resolve()),
                "sha256": sha256_file(args.plan),
            },
        }
    )
    _write_json(args.output, output)
    return 0


def _localization_digest(record: dict[str, Any]) -> dict[str, Any]:
    spin = record.get("spin_density", {})
    return {
        "radius_of_gyration_angstrom": spin.get("radius_of_gyration_angstrom"),
        "positive_spin_fraction_by_component": record.get(
            "positive_spin_fraction_by_component", {}
        ),
        "signed_spin_integral": spin.get("signed_integral"),
    }


def run_summary(args: argparse.Namespace) -> int:
    methods = _read_json(args.methods)
    settings = _settings(methods)
    production = _validate_gate(args.production_summary, args.production_gate)
    plan = _read_json(args.plan)
    if not plan.get("ready") or plan.get("scientific_status") != SCIENTIFIC_STATUS:
        raise ValueError("method-benchmark plan is not ready")
    records = [_read_json(path) for path in args.records]
    variants = [str(record["id"]) for record in settings["variants"]]
    expected = {
        (str(pair["system_id"]), int(pair["replica"]), role, variant)
        for pair in plan["selected_pairs"]
        for role in ("eda_rich", "eda_poor")
        for variant in variants
    }
    keys = [
        (
            str(record.get("system_id")),
            int(record.get("replica", -1)),
            str(record.get("seed_role")),
            str(record.get("method_variant")),
        )
        for record in records
    ]
    problems: list[str] = []
    if len(keys) != len(set(keys)):
        problems.append("duplicate method-benchmark analysis records are present")
    if set(keys) != expected:
        problems.append(
            "observed benchmark keys differ from expected: "
            f"{len(set(keys))} versus {len(expected)}"
        )
    if not records or not all(record.get("ready") for record in records):
        problems.append("one or more method-benchmark pair analyses are not ready")
    record_map = {key: record for key, record in zip(keys, records, strict=True)}
    baseline_map = {
        (str(record["system_id"]), int(record["replica"])): record
        for record in production.get("comparisons", [])
    }
    tolerance = float(settings["energy_order_tolerance_ev"])
    comparisons: list[dict[str, Any]] = []
    method_sensitive = False
    for selected in plan["selected_pairs"]:
        system = str(selected["system_id"])
        replica = int(selected["replica"])
        baseline = baseline_map.get((system, replica))
        if baseline is None:
            problems.append(f"production baseline is missing {system} r{replica}")
            continue
        baseline_delta = baseline["vertical_attachment_proxy_ev"]["rich_minus_poor"]
        baseline_delta = float(baseline_delta) if baseline_delta is not None else None
        baseline_preference = _preference(baseline_delta, tolerance)
        methods_for_pair: list[dict[str, Any]] = [
            {
                "method_variant": settings["baseline_method_id"],
                "rich_minus_poor_ev": baseline_delta,
                "preference": baseline_preference,
                "shift_from_baseline_ev": 0.0,
                "source": "accepted_stage_b2c_production",
            }
        ]
        for variant in variants:
            rich = record_map.get((system, replica, "eda_rich", variant))
            poor = record_map.get((system, replica, "eda_poor", variant))
            if rich is None or poor is None:
                continue
            rich_energy = rich["energies"]["vertical_attachment_proxy_ev"]
            poor_energy = poor["energies"]["vertical_attachment_proxy_ev"]
            delta = (
                float(rich_energy) - float(poor_energy)
                if rich_energy is not None and poor_energy is not None
                else None
            )
            shift = (
                delta - baseline_delta
                if delta is not None and baseline_delta is not None
                else None
            )
            preference = _preference(delta, tolerance)
            if preference != baseline_preference or (
                shift is not None and abs(shift) > tolerance
            ):
                method_sensitive = True
            methods_for_pair.append(
                {
                    "method_variant": variant,
                    "rich_minus_poor_ev": delta,
                    "preference": preference,
                    "shift_from_baseline_ev": shift,
                    "vertical_attachment_proxy_ev": {
                        "eda_rich": rich_energy,
                        "eda_poor": poor_energy,
                    },
                    "localization": {
                        "eda_rich": _localization_digest(rich),
                        "eda_poor": _localization_digest(poor),
                    },
                    "source": "new_method_benchmark",
                }
            )
        comparisons.append(
            {
                "system_id": system,
                "replica": replica,
                "selection_label": selected["selection_label"],
                "energy_order_tolerance_ev": tolerance,
                "methods": methods_for_pair,
            }
        )
    overall = (
        "method_sensitive_or_unresolved"
        if method_sensitive
        else "no_resolved_method_sensitivity_on_selected_pairs"
    )
    _write_json(
        args.output,
        {
            "schema_version": 1,
            "campaign": args.campaign,
            "kind": "stage_b2c_method_sensitivity_benchmark",
            "ready": not problems,
            "scientific_status": SCIENTIFIC_STATUS,
            "problems": problems,
            "overall_result": overall,
            "record_count": len(records),
            "expected_record_count": len(expected),
            "comparisons": comparisons,
            "records": records,
            "interpretation": settings["interpretation"],
            "method_limitations": [
                "representative extrema test method sensitivity; it does not resample the ensemble",
                "all calculations use fixed nuclei and are not free energies",
                "charged periodic cells remain uncorrected for finite-size artifacts",
                "PBE0 ADMM is an accuracy/cost compromise rather than a complete-basis reference",
                "method disagreement is a scientific outcome and does not make the workflow fail",
            ],
            "source": {
                "production_summary": {
                    "path": str(Path(args.production_summary).resolve()),
                    "sha256": sha256_file(args.production_summary),
                },
                "production_gate": str(Path(args.production_gate).resolve()),
                "plan": {"path": str(Path(args.plan).resolve()), "sha256": sha256_file(args.plan)},
                "methods": {
                    "path": str(Path(args.methods).resolve()),
                    "sha256": sha256_file(args.methods),
                },
            },
        },
    )
    return 0


def run_gate(args: argparse.Namespace) -> int:
    summary_path = Path(args.summary)
    summary = _read_json(summary_path)
    if not summary.get("ready"):
        raise ValueError(f"method-benchmark summary is not ready: {summary_path}")
    digest = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(f"sha256 {digest}  {summary_path.name}\n", encoding="utf-8")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan")
    plan.add_argument("--campaign", required=True)
    plan.add_argument("--production-summary", required=True)
    plan.add_argument("--production-gate", required=True)
    plan.add_argument("--methods", required=True)
    plan.add_argument("--output", required=True)
    plan.set_defaults(func=run_plan)
    render = commands.add_parser("render")
    render.add_argument("--system", required=True)
    render.add_argument("--replica", required=True, type=int)
    render.add_argument("--seed-role", required=True)
    render.add_argument("--method-variant", required=True)
    render.add_argument("--plan", required=True)
    render.add_argument("--manifest", required=True)
    render.add_argument("--methods", required=True)
    render.add_argument("--template", required=True)
    render.add_argument("--project-prefix", required=True)
    render.add_argument("--neutral-output", required=True)
    render.add_argument("--anion-output", required=True)
    render.set_defaults(func=run_render)
    analyze = commands.add_parser("analyze")
    analyze.add_argument("--campaign", required=True)
    analyze.add_argument("--system", required=True)
    analyze.add_argument("--replica", required=True, type=int)
    analyze.add_argument("--seed-role", required=True)
    analyze.add_argument("--method-variant", required=True)
    analyze.add_argument("--plan", required=True)
    analyze.add_argument("--manifest", required=True)
    analyze.add_argument("--spec", required=True)
    analyze.add_argument("--methods", required=True)
    analyze.add_argument("--systems", required=True)
    analyze.add_argument("--neutral-input", required=True)
    analyze.add_argument("--neutral-output", required=True)
    analyze.add_argument("--anion-input", required=True)
    analyze.add_argument("--anion-output", required=True)
    analyze.add_argument("--spin-cube", required=True)
    analyze.add_argument("--electron-cube", required=True)
    analyze.add_argument("--output", required=True)
    analyze.set_defaults(func=run_analyze)
    summary = commands.add_parser("summary")
    summary.add_argument("--campaign", required=True)
    summary.add_argument("--production-summary", required=True)
    summary.add_argument("--production-gate", required=True)
    summary.add_argument("--plan", required=True)
    summary.add_argument("--methods", required=True)
    summary.add_argument("--records", nargs="+", required=True)
    summary.add_argument("--output", required=True)
    summary.set_defaults(func=run_summary)
    gate = commands.add_parser("gate")
    gate.add_argument("--summary", required=True)
    gate.add_argument("--output", required=True)
    gate.set_defaults(func=run_gate)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        return int(args.func(args))
    except (KeyError, OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
