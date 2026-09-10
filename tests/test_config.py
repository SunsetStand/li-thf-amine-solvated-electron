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
        self.assertEqual(len(campaign_matrix("mixed_smoke", campaign, systems)), 1)
        self.assertEqual(len(campaign_matrix("hbond_smoke", campaign, systems)), 2)
        self.assertEqual(len(campaign_matrix("pilot", campaign, systems)), 6)
        self.assertEqual(len(campaign_matrix("hbond_pilot", campaign, systems)), 6)
        self.assertEqual(len(campaign_matrix("production", campaign, systems)), 33)

    def test_system_spec(self) -> None:
        campaign, systems, _ = load_repository_configs(ROOT)
        spec = make_system_spec("tmeda_3m", campaign, systems)
        self.assertEqual(spec.thf_count, 64)
        self.assertEqual(spec.amine_count_initial, 28)
        self.assertEqual(spec.li_electron_pairs, 1)

    def test_explicit_molar_ratio_and_neat_eda_specs(self) -> None:
        campaign, systems, _ = load_repository_configs(ROOT)
        neat = make_system_spec("pure_eda", campaign, systems)
        self.assertEqual(neat.component_counts, {"eda": 64})
        self.assertEqual(neat.thf_count, 0)
        self.assertEqual(neat.amine_count_initial, 64)
        self.assertIsNone(neat.target_concentration_m)
        self.assertEqual(neat.target_amine_mole_fraction, 1.0)

        equimolar = make_system_spec("eda_1to1", campaign, systems)
        self.assertEqual(equimolar.component_counts, {"thf": 32, "eda": 32})
        self.assertEqual(equimolar.target_amine_mole_fraction, 0.5)

    def test_repository_config_is_valid(self) -> None:
        self.assertEqual(validate_repository_configs(ROOT), [])

    def test_stage_b_mechanism_states_are_a_fixed_pair(self) -> None:
        campaign, _systems, methods = load_repository_configs(ROOT)
        state_targets = {
            state["id"]: state["li_target_valence_electrons"]
            for state in methods["stage_b_smoke"]["mechanism_states"]
        }
        self.assertEqual(state_targets, {"li0_diabatic": 3.0, "li_plus_e_diabatic": 2.0})
        self.assertEqual(campaign["campaigns"]["pilot"]["systems"], ["pure_thf", "eda_1p5m"])
        self.assertEqual(
            len(campaign["campaigns"]["pilot"]["systems"])
            * methods["stage_b"]["smoke_replicas_per_system"]
            * len(state_targets),
            4,
        )

    def test_stage_b_localization_is_explicitly_diagnostic(self) -> None:
        _campaign, systems, methods = load_repository_configs(ROOT)
        settings = methods["stage_b_localization"]
        self.assertEqual(
            settings["scientific_status"], "NUMERICAL_LOCALIZATION_DIAGNOSTIC_ONLY"
        )
        self.assertLess(settings["signed_spin_integral_min"], 1.0)
        self.assertGreater(settings["signed_spin_integral_max"], 1.0)
        for definition in [systems["thf"], *systems["amines"].values()]:
            self.assertTrue(definition["element_counts"])


if __name__ == "__main__":
    unittest.main()
