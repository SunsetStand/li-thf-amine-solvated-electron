from __future__ import annotations

import unittest

from solvelec.state import DETACHED, SOLVATED_ELECTRON, validate_state


class StateTests(unittest.TestCase):
    def test_required_states(self) -> None:
        self.assertEqual((SOLVATED_ELECTRON.charge, SOLVATED_ELECTRON.multiplicity), (0, 2))
        self.assertEqual((DETACHED.charge, DETACHED.multiplicity), (1, 1))

    def test_wrong_multiplicity_fails(self) -> None:
        with self.assertRaisesRegex(ValueError, "requires charge=0, multiplicity=2"):
            validate_state("solvated_electron", 0, 1)


if __name__ == "__main__":
    unittest.main()
