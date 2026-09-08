from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "reports" / "hbond_stage_a" / "build_report.py"


def load_report_module():
    spec = importlib.util.spec_from_file_location("hbond_stage_a_report", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class HBondStageAReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.module = load_report_module()
        cls.records = cls.module.load_records()
        cls.module.validate_records(cls.records)
        cls.metrics = cls.module.summarize(cls.records)

    def test_six_ready_solvent_only_records(self) -> None:
        self.assertEqual(len(self.records), 6)
        self.assertEqual(self.metrics["status"], "READY")
        self.assertFalse(self.metrics["scientific_scope"]["contains_li"])
        self.assertFalse(self.metrics["scientific_scope"]["contains_excess_electron"])
        self.assertTrue(self.metrics["scientific_scope"]["periodic_bulk_solvent"])

    def test_requested_molar_compositions(self) -> None:
        self.assertEqual(
            self.metrics["systems"]["pure_eda"]["composition"],
            {"thf": 0, "eda": 64, "atom_count": 768},
        )
        self.assertEqual(
            self.metrics["systems"]["eda_1to1"]["composition"],
            {"thf": 32, "eda": 32, "atom_count": 800},
        )

    def test_hydrogen_bond_results_and_cavity_boundary(self) -> None:
        pure = self.metrics["systems"]["pure_eda"]["descriptors"]
        mixed = self.metrics["systems"]["eda_1to1"]["descriptors"]
        self.assertAlmostEqual(pure["eda_eda_hydrogen_bonds"]["mean"], 79.80, places=2)
        self.assertAlmostEqual(mixed["eda_eda_hydrogen_bonds"]["mean"], 28.76, places=2)
        self.assertAlmostEqual(mixed["eda_thf_hydrogen_bonds"]["mean"], 8.04, places=2)
        self.assertGreater(pure["cavity_associated_hydrogen_bonds"]["mean"], 0.0)
        self.assertGreater(mixed["cavity_associated_hydrogen_bonds"]["mean"], 0.0)
        self.assertEqual(pure["cavity_bridging_hydrogen_bonds"]["mean"], 0.0)
        self.assertEqual(mixed["cavity_bridging_hydrogen_bonds"]["mean"], 0.0)

    def test_source_provenance(self) -> None:
        self.assertEqual(self.metrics["source_campaign"], "hbond_pilot")
        self.assertEqual(self.metrics["slurm_controller_job"], 17055)
        self.assertEqual(self.metrics["source_commit"], "969d194")


if __name__ == "__main__":
    unittest.main()
