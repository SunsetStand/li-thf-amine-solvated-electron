from __future__ import annotations

import unittest
from pathlib import Path

from solvelec.config import (
    campaign_matrix,
    format_system_id,
    load_repository_configs,
    make_system_spec,
    parse_system_id,
    validate_repository_configs,
)

ROOT = Path(__file__).resolve().parents[1]


class ConfigTests(unittest.TestCase):
    def test_system_id_round_trip(self) -> None:
        self.assertEqual(parse_system_id("eda_1p5m"), ("eda", 1.5))
        self.assertEqual(format_system_id("eda", 1.5), "eda_1p5m")
        self.assertEqual(parse_system_id("pure_thf"), (None, 0.0))

    def test_campaign_cardinality(self) -> None:
        campaign, systems, _ = load_repository_configs(ROOT)
        self.assertEqual(len(campaign_matrix("pilot", campaign, systems)), 6)
        self.assertEqual(len(campaign_matrix("production", campaign, systems)), 33)

    def test_system_spec(self) -> None:
        campaign, systems, _ = load_repository_configs(ROOT)
        spec = make_system_spec("tmeda_3m", campaign, systems)
        self.assertEqual(spec.thf_count, 64)
        self.assertEqual(spec.amine_count_initial, 28)
        self.assertEqual(spec.li_electron_pairs, 1)

    def test_repository_config_is_valid(self) -> None:
        self.assertEqual(validate_repository_configs(ROOT), [])


if __name__ == "__main__":
    unittest.main()
