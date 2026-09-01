#!/usr/bin/env python3
"""Validate solvent-only classical pilot trajectories and replica agreement."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from solvelec.composition import AVOGADRO_MOL_INV, concentration_molar, suggest_count_after_npt


def density_g_ml(total_mass_g_mol: float, volume_nm3: float) -> float:
    if total_mass_g_mol <= 0 or volume_nm3 <= 0:
        raise ValueError("mass and volume must be positive")
    volume_ml = volume_nm3 * 1.0e-21
    return total_mass_g_mol / AVOGADRO_MOL_INV / volume_ml


def replica_metrics(
    *,
    total_mass_g_mol: float,
    volumes_nm3: list[float],
    times_ps: list[float],
    amine_count: int,
    target_concentration_m: float,
    expected_duration_ns: float,
    concentration_tolerance_m: float,
    minimum_trajectory_fraction: float,
    density_half_relative_tolerance: float,
    engine_converged: bool,
) -> dict[str, Any]:
    volumes = np.asarray(volumes_nm3, dtype=float)
    times = np.asarray(times_ps, dtype=float)
    if volumes.ndim != 1 or times.ndim != 1 or len(volumes) != len(times):
        raise ValueError("trajectory volumes and times must be same-length vectors")
    if len(volumes) < 4 or not np.all(np.isfinite(volumes)) or np.any(volumes <= 0):
        raise ValueError("trajectory must contain at least four finite positive volumes")
    if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError("trajectory times must be finite and strictly increasing")
    if amine_count < 0 or target_concentration_m < 0 or expected_duration_ns <= 0:
        raise ValueError("counts and target values are invalid")

    mean_volume = float(np.mean(volumes))
    mean_density = density_g_ml(total_mass_g_mol, mean_volume)
    densities = total_mass_g_mol / AVOGADRO_MOL_INV / (volumes * 1.0e-21)
    split = len(densities) // 2
    first_half_density = float(np.mean(densities[:split]))
    second_half_density = float(np.mean(densities[split:]))
    density_half_relative_difference = abs(first_half_density - second_half_density) / float(
        np.mean(densities)
    )
    sampled_duration_ns = float((times[-1] - times[0]) / 1000.0)
    achieved_concentration = concentration_molar(amine_count, mean_volume)
    concentration_error = abs(achieved_concentration - target_concentration_m)

    checks = {
        "engine_converged": bool(engine_converged),
        "trajectory_complete": sampled_duration_ns
        >= expected_duration_ns * minimum_trajectory_fraction,
        "density_half_stable": density_half_relative_difference <= density_half_relative_tolerance,
        "concentration_within_tolerance": concentration_error <= concentration_tolerance_m,
    }
    return {
        "ready": all(checks.values()),
        "checks": checks,
        "frame_count": int(len(volumes)),
        "sampled_duration_ns": sampled_duration_ns,
        "expected_duration_ns": expected_duration_ns,
        "mean_volume_nm3": mean_volume,
        "volume_std_nm3": float(np.std(volumes, ddof=1)),
        "mean_density_g_ml": mean_density,
        "density_std_g_ml": float(np.std(densities, ddof=1)),
        "first_half_density_g_ml": first_half_density,
        "second_half_density_g_ml": second_half_density,
        "density_half_relative_difference": density_half_relative_difference,
        "achieved_concentration_m": achieved_concentration,
        "target_concentration_m": target_concentration_m,
        "concentration_error_m": concentration_error,
        "suggested_amine_count": suggest_count_after_npt(
            target_concentration_m, mean_volume, amine_count
        ),
    }


def summarize_records(
    records: list[dict[str, Any]], replica_density_relative_tolerance: float
) -> dict[str, Any]:
    if not records or replica_density_relative_tolerance <= 0:
        raise ValueError("records and a positive replica tolerance are required")
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[str(record["system_id"])].append(record)

    systems: dict[str, Any] = {}
    for system_id, group in sorted(grouped.items()):
        available_densities = [
            float(row["metrics"]["mean_density_g_ml"])
            for row in group
            if "mean_density_g_ml" in row.get("metrics", {})
        ]
        densities = np.asarray(available_densities)
        if len(densities):
            mean_density: float | None = float(np.mean(densities))
            replica_relative_span: float | None = (
                float((np.max(densities) - np.min(densities)) / mean_density)
                if len(densities) > 1
                else 0.0
            )
        else:
            mean_density = None
            replica_relative_span = None
        checks = {
            "all_replicas_ready": all(bool(row["metrics"]["ready"]) for row in group),
            "replica_density_consistent": replica_relative_span is not None
            and len(densities) == len(group)
            and replica_relative_span <= replica_density_relative_tolerance,
        }
        systems[system_id] = {
            "ready": all(checks.values()),
            "checks": checks,
            "replicas": sorted(int(row["replica"]) for row in group),
            "mean_density_g_ml": mean_density,
            "replica_density_relative_span": replica_relative_span,
        }
    return {
        "schema_version": 1,
        "ready": all(record["ready"] for record in systems.values()),
        "replica_density_relative_tolerance": replica_density_relative_tolerance,
        "systems": systems,
    }


def _read_trajectory(tpr: Path, trajectory: Path) -> tuple[float, list[float], list[float]]:
    import MDAnalysis as mda

    universe = mda.Universe(str(tpr), str(trajectory))
    total_mass = float(universe.atoms.total_mass())
    volumes_nm3: list[float] = []
    times_ps: list[float] = []
    for frame in universe.trajectory:
        volumes_nm3.append(float(frame.volume) / 1000.0)
        times_ps.append(float(frame.time))
    return total_mass, volumes_nm3, times_ps


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run_replica(args: argparse.Namespace) -> int:
    spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
    methods = json.loads(Path(args.methods).read_text(encoding="utf-8"))
    engine = json.loads(Path(args.engine_validation).read_text(encoding="utf-8"))
    try:
        mass, volumes, times = _read_trajectory(Path(args.tpr), Path(args.trajectory))
        metrics = replica_metrics(
            total_mass_g_mol=mass,
            volumes_nm3=volumes,
            times_ps=times,
            amine_count=int(spec["amine_count_initial"]),
            target_concentration_m=float(spec["target_concentration_m"]),
            expected_duration_ns=float(methods["classical_md"]["production_ns"]),
            concentration_tolerance_m=float(args.concentration_tolerance_m),
            minimum_trajectory_fraction=float(
                methods["classical_validation"]["minimum_trajectory_fraction"]
            ),
            density_half_relative_tolerance=float(
                methods["classical_validation"]["density_half_relative_tolerance"]
            ),
            engine_converged=bool(engine.get("converged")),
        )
        result = {
            "schema_version": 1,
            "system_id": spec["system_id"],
            "replica": int(spec["replica"]),
            "solvent_only": True,
            "metrics": metrics,
            "inputs": {
                "spec": str(Path(args.spec).resolve()),
                "tpr": str(Path(args.tpr).resolve()),
                "trajectory": str(Path(args.trajectory).resolve()),
                "engine_validation": str(Path(args.engine_validation).resolve()),
            },
        }
    except Exception as exc:
        # Preserve a machine-readable failed gate for corrupt or unsupported
        # trajectories; the later gate rule will stop the campaign.
        result = {
            "schema_version": 1,
            "system_id": spec.get("system_id"),
            "replica": spec.get("replica"),
            "solvent_only": True,
            "metrics": {"ready": False, "error": str(exc)},
        }
    _write_json(Path(args.output), result)
    return 0


def run_summary(args: argparse.Namespace) -> int:
    methods = json.loads(Path(args.methods).read_text(encoding="utf-8"))
    records = [json.loads(Path(path).read_text(encoding="utf-8")) for path in args.inputs]
    result = summarize_records(
        records,
        float(methods["classical_validation"]["replica_density_relative_tolerance"]),
    )
    result["campaign"] = args.campaign
    _write_json(Path(args.output), result)
    return 0


def run_gate(args: argparse.Namespace) -> int:
    summary_path = Path(args.summary)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    output = Path(args.output)
    output.unlink(missing_ok=True)
    if not summary.get("ready"):
        print(f"ERROR: classical pilot validation failed; inspect {summary_path}", file=sys.stderr)
        return 4
    digest = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(f"summary_sha256={digest}\n", encoding="utf-8")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    replica = subparsers.add_parser("replica")
    replica.add_argument("--spec", required=True)
    replica.add_argument("--methods", required=True)
    replica.add_argument("--tpr", required=True)
    replica.add_argument("--trajectory", required=True)
    replica.add_argument("--engine-validation", required=True)
    replica.add_argument("--concentration-tolerance-m", type=float, required=True)
    replica.add_argument("--output", required=True)
    replica.set_defaults(func=run_replica)

    summary = subparsers.add_parser("summary")
    summary.add_argument("--campaign", required=True)
    summary.add_argument("--methods", required=True)
    summary.add_argument("--output", required=True)
    summary.add_argument("inputs", nargs="+")
    summary.set_defaults(func=run_summary)

    gate = subparsers.add_parser("gate")
    gate.add_argument("--summary", required=True)
    gate.add_argument("--output", required=True)
    gate.set_defaults(func=run_gate)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
