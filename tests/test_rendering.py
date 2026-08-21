from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from solvelec.config import load_repository_configs
from solvelec.rendering import render_cp2k, render_packmol

ROOT = Path(__file__).resolve().parents[1]


class RenderingTests(unittest.TestCase):
    def test_cp2k_state_and_constraint(self) -> None:
        _, _, methods = load_repository_configs(ROOT)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "test.inp"
            render_cp2k(
                ROOT / "workflow" / "templates" / "cp2k" / "pbe0_cdft.inp.tpl",
                output,
                "solvated_electron",
                "smoke",
                "frame.xyz",
                "cell.inc",
                methods["cp2k"],
                li_atom_index=1,
                constrained=True,
            )
            text = output.read_text(encoding="utf-8")
            self.assertIn("CHARGE 0", text)
            self.assertIn("MULTIPLICITY 2", text)
            self.assertIn("&CDFT", text)
            self.assertIn("ATOMS 1", text)
            self.assertIn("&BECKE_CONSTRAINT", text)
            self.assertIn("ADMM_TYPE ADMMS", text)

    def test_detached_cannot_use_li_constraint(self) -> None:
        _, _, methods = load_repository_configs(ROOT)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                render_cp2k(
                    ROOT / "workflow" / "templates" / "cp2k" / "pbe0_cdft.inp.tpl",
                    Path(directory) / "bad.inp",
                    "detached",
                    "bad",
                    "frame.xyz",
                    "cell.inc",
                    methods["cp2k"],
                    li_atom_index=1,
                    constrained=True,
                )

    def test_packmol_counts(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "packmol.inp"
            render_packmol(output, "box.pdb", 22.0, "thf.pdb", 64, "eda.pdb", 9, seed=2)
            text = output.read_text(encoding="utf-8")
            self.assertIn("number 64", text)
            self.assertIn("number 9", text)
            self.assertIn("seed 2", text)


if __name__ == "__main__":
    unittest.main()
