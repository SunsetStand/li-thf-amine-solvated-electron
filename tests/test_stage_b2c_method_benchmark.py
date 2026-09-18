from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from solvelec.rendering import render_stage_b2c_method_benchmark_cp2k

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "stage_b2c_method_benchmark.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("stage_b2c_method_benchmark", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _write_gate(summary: Path, gate: Path) -> None:
    digest = hashlib.sha256(summary.read_bytes()).hexdigest()
    gate.write_text(f"sha256 {digest}  {summary.name}\n", encoding="utf-8")


def _comparison(system: str, replica: int, delta: float) -> dict:
    return {
        "system_id": system,
        "replica": replica,
        "ready": True,
        "vertical_attachment_proxy_ev": {
            "eda_rich": 1.0 + delta,
            "eda_poor": 1.0,
            "rich_minus_poor": delta,
        },
    }


class StageB2CMethodBenchmarkTests(unittest.TestCase):
    def test_tight_pbe_and_pbe0_inputs_are_method_specific(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            coordinates = root / "coordinates.xyz"
            cell = root / "cell.inc"
            coordinates.write_text("2\ntest\nH 1 1 1\nGh 5 5 5\n", encoding="utf-8")
            cell.write_text(
                "&CELL\n  A 10 0 0\n  B 0 10 0\n  C 0 0 10\n"
                "  PERIODIC XYZ\n&END CELL\n",
                encoding="utf-8",
            )
            methods = json.loads((ROOT / "configs" / "methods.yaml").read_text())
            base = methods["stage_b2c_preferential_smoke"]
            states = {record["id"]: record for record in base["states"]}
            for variant in methods["stage_b2c_method_benchmark"]["variants"]:
                method = {
                    **base,
                    **variant,
                    "scientific_status": "NUMERICAL_METHOD_SENSITIVITY_BENCHMARK_ONLY",
                }
                for state_id in ("neutral", "anion"):
                    output = root / variant["id"] / state_id / "cp2k.inp"
                    render_stage_b2c_method_benchmark_cp2k(
                        ROOT
                        / "workflow"
                        / "templates"
                        / "cp2k"
                        / "stage_b2c_method_benchmark.inp.tpl",
                        output,
                        project="test",
                        coordinates_path=coordinates,
                        cell_path=cell,
                        method=method,
                        state=states[state_id],
                    )
                    text = output.read_text(encoding="utf-8")
                    self.assertIn("RUN_TYPE ENERGY", text)
                    self.assertIn("TZV2P-MOLOPT-GTH", text)
                    self.assertNotIn("&KIND LI", text.upper())
                    self.assertNotIn("&CDFT", text.upper())
                    self.assertEqual("&E_DENSITY_CUBE" in text, state_id == "anion")
                    self.assertEqual("&HF" in text, variant["id"] == "pbe0_admm")
                    self.assertEqual(
                        "AUXILIARY_DENSITY_MATRIX_METHOD" in text,
                        variant["id"] == "pbe0_admm",
                    )

    def test_plan_selects_and_locks_production_extrema(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            summary = root / "production.summary.json"
            gate = root / "production.done"
            _write_json(
                summary,
                {
                    "ready": True,
                    "scientific_status": (
                        "NUMERICAL_PREFERENTIAL_SOLVATION_ENSEMBLE_SCREEN_ONLY"
                    ),
                    "comparisons": [
                        _comparison("eda_1p5m", 1, 0.01),
                        _comparison("eda_1p5m", 2, 0.04),
                        _comparison("eda_1p5m", 3, -0.03),
                        _comparison("eda_3m", 1, 0.07),
                        _comparison("eda_3m", 2, 0.02),
                        _comparison("eda_3m", 3, -0.04),
                    ],
                },
            )
            _write_gate(summary, gate)
            output = root / "plan.json"
            self.assertEqual(
                module.run_plan(
                    SimpleNamespace(
                        campaign="pilot",
                        production_summary=str(summary),
                        production_gate=str(gate),
                        methods=str(ROOT / "configs" / "methods.yaml"),
                        output=str(output),
                    )
                ),
                0,
            )
            plan = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(plan["ready"])
            self.assertEqual(plan["expected_cp2k_single_points"], 32)
            self.assertEqual(
                [(record["system_id"], record["replica"]) for record in plan["selected_pairs"]],
                [("eda_1p5m", 2), ("eda_1p5m", 3), ("eda_3m", 1), ("eda_3m", 3)],
            )

    def test_summary_records_method_disagreement_as_scientific_result(self) -> None:
        module = _load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            production = root / "production.summary.json"
            gate = root / "production.done"
            comparisons = [
                _comparison("eda_1p5m", 2, 0.04),
                _comparison("eda_1p5m", 3, -0.03),
                _comparison("eda_3m", 1, 0.07),
                _comparison("eda_3m", 3, -0.04),
            ]
            _write_json(
                production,
                {
                    "ready": True,
                    "scientific_status": (
                        "NUMERICAL_PREFERENTIAL_SOLVATION_ENSEMBLE_SCREEN_ONLY"
                    ),
                    "comparisons": comparisons,
                },
            )
            _write_gate(production, gate)
            plan = root / "plan.json"
            selected = []
            for comparison in comparisons:
                selected.append(
                    {
                        "system_id": comparison["system_id"],
                        "replica": comparison["replica"],
                        "selection_label": "extremum",
                        "baseline_comparison": comparison,
                    }
                )
            _write_json(
                plan,
                {
                    "ready": True,
                    "scientific_status": "NUMERICAL_METHOD_SENSITIVITY_BENCHMARK_ONLY",
                    "selected_pairs": selected,
                },
            )
            record_paths: list[str] = []
            for pair in selected:
                for variant in ("pbe_converged", "pbe0_admm"):
                    delta = 0.03
                    if pair["system_id"] == "eda_1p5m" and pair["replica"] == 2:
                        delta = -0.08 if variant == "pbe0_admm" else 0.04
                    for role, energy in (("eda_rich", 1.0 + delta), ("eda_poor", 1.0)):
                        path = root / (
                            f"{pair['system_id']}-r{pair['replica']}-{role}-{variant}.json"
                        )
                        _write_json(
                            path,
                            {
                                "ready": True,
                                "system_id": pair["system_id"],
                                "replica": pair["replica"],
                                "seed_role": role,
                                "method_variant": variant,
                                "energies": {"vertical_attachment_proxy_ev": energy},
                                "positive_spin_fraction_by_component": {"eda": 0.2},
                                "spin_density": {
                                    "signed_integral": 1.0,
                                    "radius_of_gyration_angstrom": 4.0,
                                },
                            },
                        )
                        record_paths.append(str(path))
            output = root / "benchmark.summary.json"
            self.assertEqual(
                module.run_summary(
                    SimpleNamespace(
                        campaign="pilot",
                        production_summary=str(production),
                        production_gate=str(gate),
                        plan=str(plan),
                        methods=str(ROOT / "configs" / "methods.yaml"),
                        records=record_paths,
                        output=str(output),
                    )
                ),
                0,
            )
            summary = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(summary["ready"])
            self.assertEqual(summary["record_count"], 16)
            self.assertEqual(summary["overall_result"], "method_sensitive_or_unresolved")
            self.assertEqual(summary["problems"], [])

    def test_workflow_is_32_cp2k_single_points_without_gromacs(self) -> None:
        methods = json.loads((ROOT / "configs" / "methods.yaml").read_text())
        settings = methods["stage_b2c_method_benchmark"]
        pair_count = sum(len(values) for values in settings["representative_replicas"].values())
        self.assertEqual(pair_count * 2 * 2 * len(settings["variants"]), 32)
        rules = (ROOT / "workflow" / "rules" / "62_stage_b2c_method_benchmark.smk").read_text()
        self.assertEqual(rules.count("bash {STAGE_RUNNER:q} cdft --"), 2)
        self.assertNotIn("gmx", rules.lower())
        self.assertNotIn("&CDFT", rules.upper())
        self.assertIn("stage_b2c_method_benchmark", (ROOT / "run.sh").read_text())


if __name__ == "__main__":
    unittest.main()
