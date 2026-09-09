from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from solvelec.parsers import HARTREE_TO_EV

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "prepare_stage_b.py"


def load_script():
    spec = importlib.util.spec_from_file_location("prepare_stage_b", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cp2k_output(target: float, current: float, energy: float) -> str:
    return f"""
 CDFT SCF iter = 4 RMS gradient = 0.01 energy = {energy:.12f}
 Target value of constraint  : {target:.12f}
 Current value of constraint : {current:.12f}
 Deviation from target       : {current - target:.12f}
 Strength of constraint      : 1.250000000000
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: {energy:.12f}
 PROGRAM ENDED AT 2026-09-08 00:00:00
 """


class StageBMechanismSmokeTests(unittest.TestCase):
    def test_paired_summary_reports_fixed_geometry_energy_difference(self) -> None:
        module = load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "ready": True,
                        "system_id": "eda_1p5m",
                        "amine": "eda",
                        "replica": 1,
                        "candidates": [{"candidate_id": "separated", "ready": True}],
                    }
                ),
                encoding="utf-8",
            )
            inputs = [root / "li0.inp", root / "li_plus_e.inp"]
            outputs = [root / "li0.out", root / "li_plus_e.out"]
            for path in inputs:
                path.write_text("synthetic input\n", encoding="utf-8")
            outputs[0].write_text(cp2k_output(3.0, 3.01, -100.0), encoding="utf-8")
            outputs[1].write_text(cp2k_output(2.0, 2.02, -99.9), encoding="utf-8")
            summary = root / "summary.json"
            args = SimpleNamespace(
                campaign="pilot",
                candidate="separated",
                methods=str(ROOT / "configs" / "methods.yaml"),
                output=str(summary),
                outputs=[str(path) for path in outputs],
                cp2k_inputs=[str(path) for path in inputs],
                manifests=[str(manifest), str(manifest)],
                states=["li0_diabatic", "li_plus_e_diabatic"],
                targets=[3.0, 2.0],
            )
            self.assertEqual(module.run_mechanism_summary(args), 0)
            result = json.loads(summary.read_text(encoding="utf-8"))
            self.assertTrue(result["ready"])
            self.assertIn("NOT_A_STABILITY_OR_LOCALIZATION_RESULT", result["scientific_status"])
            self.assertEqual(len(result["records"]), 1)
            record = result["records"][0]
            self.assertTrue(record["complete_pair"])
            self.assertAlmostEqual(record["delta_energy_hartree_li_plus_e_minus_li0"], 0.1)
            self.assertAlmostEqual(record["fixed_geometry_diabatic_gap_ev"], 0.1 * HARTREE_TO_EV)
            self.assertEqual(
                {state["li_target_valence_electrons"] for state in record["states"]},
                {2.0, 3.0},
            )

    def test_paired_summary_fails_ready_gate_when_one_constraint_misses(self) -> None:
        module = load_script()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "ready": True,
                        "system_id": "pure_thf",
                        "amine": None,
                        "replica": 1,
                        "candidates": [{"candidate_id": "separated", "ready": True}],
                    }
                ),
                encoding="utf-8",
            )
            inputs = [root / "li0.inp", root / "li_plus_e.inp"]
            outputs = [root / "li0.out", root / "li_plus_e.out"]
            for path in inputs:
                path.write_text("synthetic input\n", encoding="utf-8")
            outputs[0].write_text(cp2k_output(3.0, 3.01, -100.0), encoding="utf-8")
            outputs[1].write_text(cp2k_output(2.0, 2.08, -99.9), encoding="utf-8")
            summary = root / "summary.json"
            args = SimpleNamespace(
                campaign="pilot",
                candidate="separated",
                methods=str(ROOT / "configs" / "methods.yaml"),
                output=str(summary),
                outputs=[str(path) for path in outputs],
                cp2k_inputs=[str(path) for path in inputs],
                manifests=[str(manifest), str(manifest)],
                states=["li0_diabatic", "li_plus_e_diabatic"],
                targets=[3.0, 2.0],
            )
            module.run_mechanism_summary(args)
            result = json.loads(summary.read_text(encoding="utf-8"))
            self.assertFalse(result["ready"])
            self.assertFalse(result["records"][0]["ready"])
            failed_state = next(
                state
                for state in result["records"][0]["states"]
                if state["state_id"] == "li_plus_e_diabatic"
            )
            self.assertIn("exceeds", failed_state["problems"][0])


if __name__ == "__main__":
    unittest.main()
