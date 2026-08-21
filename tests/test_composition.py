from __future__ import annotations

import math
import unittest

from solvelec.composition import (
    concentration_molar,
    enrichment_factor,
    initial_amine_count,
    log_odds_enrichment,
    mole_fraction,
    suggest_count_after_npt,
)


class CompositionTests(unittest.TestCase):
    def test_concentration_round_trip(self) -> None:
        concentration = concentration_molar(8, 8.854e0)
        self.assertAlmostEqual(concentration, 1.5, delta=0.01)

    def test_seed_counts_match_campaign_scale(self) -> None:
        thf_volume = 0.08120
        expected = {
            (1.5, 0.0669): 9,
            (3.0, 0.0669): 20,
            (1.5, 0.0850): 9,
            (3.0, 0.0850): 21,
            (1.5, 0.1080): 9,
            (3.0, 0.1080): 23,
            (1.5, 0.1500): 10,
            (3.0, 0.1500): 28,
        }
        for (target, additive_volume), count in expected.items():
            with self.subTest(target=target, additive_volume=additive_volume):
                self.assertEqual(
                    initial_amine_count(target, 64, thf_volume, additive_volume), count
                )

    def test_npt_count_update_is_bounded(self) -> None:
        self.assertEqual(suggest_count_after_npt(3.0, 10.0, 10, max_change=3), 13)

    def test_preferential_solvation_metrics(self) -> None:
        self.assertAlmostEqual(mole_fraction(16, 64), 0.2)
        self.assertAlmostEqual(enrichment_factor(0.4, 0.2), 2.0)
        self.assertAlmostEqual(log_odds_enrichment(0.5, 0.2), math.log(4.0))

    def test_invalid_values_fail(self) -> None:
        with self.assertRaises(ValueError):
            concentration_molar(1, 0)
        with self.assertRaises(ValueError):
            initial_amine_count(10.0, 64, 0.0812, 0.15)


if __name__ == "__main__":
    unittest.main()
