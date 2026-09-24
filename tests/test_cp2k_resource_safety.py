from __future__ import annotations

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RULE_PATHS = (
    ROOT / "workflow" / "rules" / "55_stage_b.smk",
    ROOT / "workflow" / "rules" / "56_stage_b_mechanism_smoke.smk",
    ROOT / "workflow" / "rules" / "59_stage_b2_intrinsic.smk",
    ROOT / "workflow" / "rules" / "61_stage_b2c_preferential.smk",
    ROOT / "workflow" / "rules" / "62_stage_b2c_method_benchmark.smk",
    ROOT / "workflow" / "rules" / "63_stage_b3_amine_series.smk",
)
PROFILE_PATHS = (
    ROOT / "configs" / "profiles" / "slurm" / "config.v9+.yaml",
    ROOT / "configs" / "profiles" / "tmc-amd" / "config.v9+.yaml",
)


class Cp2kResourceSafetyTests(unittest.TestCase):
    def test_cp2k_rules_fit_the_twelve_cpu_child_budget(self) -> None:
        combined = "\n".join(path.read_text(encoding="utf-8") for path in RULE_PATHS)
        self.assertNotIn("tasks=32", combined)
        self.assertNotIn("tasks=8", combined)
        self.assertEqual(combined.count("tasks=12,"), 9)
        self.assertEqual(combined.count("cp2k_slots=1,"), 9)
        self.assertEqual(combined.count("cpu_slots=12,"), 9)

    def test_slurm_profiles_cap_jobs_and_cp2k_concurrency(self) -> None:
        for path in PROFILE_PATHS:
            with self.subTest(profile=path.parent.name):
                profile = path.read_text(encoding="utf-8")
                self.assertIn("jobs: 8", profile)
                self.assertIn("resources:\n  cp2k_slots: 1", profile)
                self.assertIn("  cpu_slots: 12", profile)
                self.assertIn("  - cpu_slots=4", profile)
                self.assertNotIn("tasks: 32", profile)
                self.assertNotIn("tasks: 8", profile)
                self.assertEqual(profile.count("tasks: 12"), 11)


if __name__ == "__main__":
    unittest.main()
