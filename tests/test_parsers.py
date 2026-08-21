from __future__ import annotations

import unittest

from solvelec.parsers import parse_cp2k_text, parse_orca_text, vertical_detachment_energy_ev


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

    def test_orca_success(self) -> None:
        text = "FINAL SINGLE POINT ENERGY -10.25\nORCA TERMINATED NORMALLY\n"
        self.assertTrue(parse_orca_text(text).converged)

    def test_vde_sign(self) -> None:
        self.assertAlmostEqual(vertical_detachment_energy_ev(-100.0, -99.9), 2.7211386246)


if __name__ == "__main__":
    unittest.main()
