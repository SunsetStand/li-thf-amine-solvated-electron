from __future__ import annotations

import unittest

from solvelec.parsers import (
    evaluate_cp2k_cdft_constraint,
    parse_cp2k_cdft_iterations,
    parse_cp2k_text,
    parse_gromacs_text,
    parse_orca_text,
    parse_packmol_text,
    vertical_detachment_energy_ev,
)


class ParserTests(unittest.TestCase):
    def test_cp2k_success(self) -> None:
        text = """
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -123.456789
 PROGRAM ENDED AT 2026-08-21 00:00:00
 """
        result = parse_cp2k_text(text)
        self.assertTrue(result.converged)
        self.assertAlmostEqual(result.energy_hartree or 0, -123.456789)

    def test_cp2k_failure_is_not_success(self) -> None:
        text = "SCF run NOT converged\nPROGRAM ENDED AT now\n"
        result = parse_cp2k_text(text)
        self.assertFalse(result.converged)
        self.assertIsNone(result.energy_hartree)

    def test_cp2k_cdft_gate_uses_final_population_and_configured_tolerance(self) -> None:
        text = """
 CDFT SCF iter =     3 RMS gradient =   0.95E+00 energy = -3117.7730892854
 Target value of constraint  : 2.000000000000
 Current value of constraint : 1.047863017043
 Deviation from target       : -9.521E-01
 Strength of constraint      : 3.087628454047
 CDFT SCF iter =     4 RMS gradient =   0.45E-01 energy = -3117.3215608767
 Target value of constraint  : 2.000000000000
 Current value of constraint : 2.045021000000
 Deviation from target       : 4.502E-02
 Strength of constraint      : 1.258500000000
 """
        iterations = parse_cp2k_cdft_iterations(text)
        self.assertEqual(len(iterations), 2)
        self.assertEqual(iterations[-1].iteration, 4)
        relaxed = evaluate_cp2k_cdft_constraint(
            text, expected_target_electrons=2.0, tolerance_electrons=0.05
        )
        strict = evaluate_cp2k_cdft_constraint(
            text, expected_target_electrons=2.0, tolerance_electrons=0.02
        )
        self.assertTrue(relaxed.converged)
        self.assertAlmostEqual(relaxed.deviation_electrons or 0.0, 0.045021)
        self.assertFalse(strict.converged)
        self.assertIn("exceeds", strict.problems[0])

    def test_cp2k_cdft_gate_rejects_missing_iteration_and_wrong_target(self) -> None:
        missing = evaluate_cp2k_cdft_constraint(
            "PROGRAM ENDED AT now", expected_target_electrons=2.0, tolerance_electrons=0.05
        )
        self.assertFalse(missing.converged)
        self.assertIn("no completed", missing.problems[0])
        wrong_target_text = """
 CDFT SCF iter = 1 RMS gradient = 0.01 energy = -1.0
 Target value of constraint: 1.0
 Current value of constraint: 1.0
 Deviation from target: 0.0
 Strength of constraint: 0.5
 """
        wrong_target = evaluate_cp2k_cdft_constraint(
            wrong_target_text, expected_target_electrons=2.0, tolerance_electrons=0.05
        )
        self.assertFalse(wrong_target.converged)
        self.assertIn("differs", wrong_target.problems[0])

    def test_orca_success(self) -> None:
        text = "FINAL SINGLE POINT ENERGY -10.25\nORCA TERMINATED NORMALLY\n"
        self.assertTrue(parse_orca_text(text).converged)

    def test_vde_sign(self) -> None:
        self.assertAlmostEqual(vertical_detachment_energy_ev(-100.0, -99.9), 2.7211386246)

    def test_packmol_requires_success_without_error(self) -> None:
        success = parse_packmol_text("Success!\nFinal objective function value: 0.1")
        self.assertTrue(success.converged)
        self.assertFalse(parse_packmol_text("ERROR: could not open molecule.pdb").converged)

    def test_gromacs_requires_finished_marker_without_fatal_error(self) -> None:
        self.assertTrue(parse_gromacs_text("Finished mdrun on rank 0").converged)
        self.assertFalse(parse_gromacs_text("Fatal error: bad topology").converged)


if __name__ == "__main__":
    unittest.main()
