from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "reports" / "stage_a" / "build_report.py"


def load_report_module():
    spec = importlib.util.spec_from_file_location("stage_a_report", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class StageAReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.module = load_report_module()
        cls.records = cls.module.load_records()
        cls.module.validate_records(cls.records)
        cls.metrics = cls.module.summarize(cls.records)

    def test_committed_snapshot_bank_is_ready_and_solvent_only(self) -> None:
        self.assertEqual(len(self.records), 6)
        self.assertEqual(self.metrics["status"], "READY")
        self.assertFalse(self.metrics["scientific_scope"]["contains_li"])
        self.assertFalse(self.metrics["scientific_scope"]["contains_excess_electron"])

    def test_mixed_composition_and_concentration(self) -> None:
        mixed = self.metrics["systems"]["eda_1p5m"]
        self.assertEqual(mixed["thf_count"], 64)
        self.assertEqual(mixed["eda_count"], 9)
        self.assertAlmostEqual(mixed["thf_to_eda_ratio"], 64 / 9)
        self.assertAlmostEqual(mixed["eda_mole_fraction"], 9 / 73)
        self.assertAlmostEqual(mixed["achieved_eda_concentration_m"], 1.5102, places=4)

    def test_reported_ensemble_density_and_void_metrics(self) -> None:
        pure = self.metrics["systems"]["pure_thf"]
        mixed = self.metrics["systems"]["eda_1p5m"]
        self.assertAlmostEqual(pure["density_g_ml"]["mean"], 0.87697924, places=8)
        self.assertAlmostEqual(mixed["density_g_ml"]["mean"], 0.86569301, places=8)
        self.assertAlmostEqual(pure["void_radius_angstrom"]["mean"], 1.89385220, places=8)
        self.assertAlmostEqual(mixed["void_radius_angstrom"]["mean"], 1.97047491, places=8)


if __name__ == "__main__":
    unittest.main()
