from __future__ import annotations

import unittest

from solvelec.classification import LocalizationMetrics, classify_localization

THRESHOLDS = {
    "electron_count_min": 0.75,
    "electron_count_max": 1.25,
    "li_spin_collapse": 0.50,
    "atomic_spin_localized": 0.50,
    "cavity_interstitial_fraction": 0.55,
}


class ClassificationTests(unittest.TestCase):
    def test_cavity(self) -> None:
        result = classify_localization(LocalizationMetrics(1.0, 0.05, 0.10, 0.80), THRESHOLDS)
        self.assertEqual(result.label, "cavity_electron")

    def test_li_collapse_has_priority(self) -> None:
        result = classify_localization(LocalizationMetrics(1.0, 0.80, 0.80, 0.80), THRESHOLDS)
        self.assertEqual(result.label, "li_atomic_or_contact")

    def test_molecular_anion(self) -> None:
        result = classify_localization(LocalizationMetrics(1.0, 0.05, 0.70, 0.20), THRESHOLDS)
        self.assertEqual(result.label, "molecular_anion")

    def test_bad_spin_integral(self) -> None:
        result = classify_localization(LocalizationMetrics(0.3, 0.0, 0.0, 1.0), THRESHOLDS)
        self.assertEqual(result.label, "invalid_spin_integral")


if __name__ == "__main__":
    unittest.main()
