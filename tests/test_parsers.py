from __future__ import annotations

import unittest

from solvelec.parsers import (
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
