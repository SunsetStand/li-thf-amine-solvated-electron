from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path
from types import ModuleType

import numpy as np

from solvelec.trajectory import (
    autocorrelation_summary,
    cell_matrix,
    largest_void_proxy,
    pair_distances,
    select_representative_indices,
)

ROOT = Path(__file__).resolve().parents[1]


def load_analysis_script() -> ModuleType:
    path = ROOT / "workflow" / "scripts" / "analyze_classical_ensemble.py"
    spec = importlib.util.spec_from_file_location("analyze_classical_ensemble", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TrajectoryAnalysisTests(unittest.TestCase):
    def test_cell_and_pair_distance_cross_periodic_boundary(self) -> None:
        cell = cell_matrix([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
        np.testing.assert_allclose(cell, np.diag([10.0, 10.0, 10.0]), atol=1.0e-12)
        distances = pair_distances(
            np.asarray([[0.5, 5.0, 5.0]]),
            np.asarray([[9.5, 5.0, 5.0]]),
            cell,
        )
        np.testing.assert_allclose(distances, [1.0], atol=1.0e-12)

    def test_same_residue_pairs_can_be_excluded(self) -> None:
        cell = np.diag([10.0, 10.0, 10.0])
        positions = np.asarray([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [5.0, 1.0, 1.0]])
        residues = np.asarray([1, 1, 2])
        distances = pair_distances(
            positions,
            positions,
            cell,
            same_group=True,
            first_residue_ids=residues,
            second_residue_ids=residues,
            exclude_same_residue=True,
        )
        np.testing.assert_allclose(np.sort(distances), [3.0, 4.0])

    def test_void_proxy_finds_open_region(self) -> None:
        result = largest_void_proxy(
            np.asarray([[5.0, 5.0, 5.0]]),
            np.asarray([1.0]),
            np.diag([10.0, 10.0, 10.0]),
            points_per_axis=6,
            refinement_levels=2,
        )
        self.assertGreater(result["radius_angstrom"], 7.0)
        self.assertEqual(len(result["fractional"]), 3)

    def test_autocorrelation_reports_effective_samples(self) -> None:
        times = np.arange(20, dtype=float) * 100.0
        constant = autocorrelation_summary(times, np.ones(20))
        self.assertEqual(constant["effective_sample_size"], 20.0)
        correlated = autocorrelation_summary(times, np.linspace(0.0, 1.0, 20))
        self.assertGreater(correlated["statistical_inefficiency"], 1.0)
        self.assertLess(correlated["effective_sample_size"], 20.0)

    def test_snapshot_selection_is_deterministic_and_equilibrated(self) -> None:
        records = [
            {"time_ps": float(index * 1000), "density": 0.8 + index * 0.01, "void": 2.0}
            for index in range(10)
        ]
        first = select_representative_indices(
            records,
            ["density", "void"],
            count=2,
            minimum_time_ps=4000.0,
            minimum_separation_ps=2000.0,
        )
        second = select_representative_indices(
            records,
            ["density", "void"],
            count=2,
            minimum_time_ps=4000.0,
            minimum_separation_ps=2000.0,
        )
        self.assertEqual(first, second)
        self.assertTrue(all(records[index]["time_ps"] >= 4000.0 for index in first))
        self.assertGreaterEqual(
            abs(records[first[0]]["time_ps"] - records[first[1]]["time_ps"]), 2000.0
        )

    def test_failed_analysis_record_survives_summary(self) -> None:
        module = load_analysis_script()
        records = [
            {
                "system_id": "pure_thf",
                "replica": 1,
                "ready": False,
                "metrics": {"ready": False, "error": "trajectory unreadable"},
            }
        ]
        summary = module.summarize_records(records, "analysis")
        self.assertFalse(summary["ready"])
        self.assertIsNone(summary["systems"]["pure_thf"]["minimum_effective_samples"])


if __name__ == "__main__":
    unittest.main()
