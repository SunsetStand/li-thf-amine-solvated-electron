from __future__ import annotations

import json
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from solvelec.cli import main
from solvelec.config import load_repository_configs
from solvelec.rendering import render_cp2k, render_gromacs_mdp, render_packmol

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
            self.assertIn("TARGET 2.0", text)
            self.assertIn("POTENTIAL GTH-PBE-q3", text)
            self.assertIn("&BECKE_CONSTRAINT", text)
            self.assertIn("ADMM_TYPE ADMMS", text)
            self.assertIn("&E_DENSITY_CUBE", text)
            self.assertNotIn("&SPIN_DENSITY_CUBE", text)

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

    def test_packmol_cli_accepts_explicit_production_paths(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "packmol.inp"
            with (
                patch("solvelec.cli.repository_root", return_value=ROOT),
                redirect_stdout(StringIO()),
            ):
                result = main(
                    [
                        "render-packmol",
                        "--system",
                        "pure_thf",
                        "--replica",
                        "1",
                        "--output",
                        str(output),
                        "--output-pdb",
                        "/storage/packed.pdb",
                        "--thf-structure",
                        "/storage/thf.pdb",
                    ]
                )

            self.assertEqual(result, 0)
            text = output.read_text(encoding="utf-8")
            self.assertIn("output /storage/packed.pdb", text)
            self.assertIn("structure /storage/thf.pdb", text)

    def test_smoke_mdp_has_reproducible_thermodynamic_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            template = directory / "nvt.mdp.tpl"
            output = directory / "nvt.mdp"
            template.write_text(
                "ref-t = $temperature_k\nref-p = $pressure_bar\ngen-seed = $seed\n",
                encoding="utf-8",
            )
            render_gromacs_mdp(template, output, 298.15, 1.0, 2026001)
            text = output.read_text(encoding="utf-8")
            self.assertIn("ref-t = 298.15", text)
            self.assertIn("ref-p = 1.000000", text)
            self.assertIn("gen-seed = 2026001", text)

    def test_pilot_mdp_uses_configured_duration_and_stable_system_seed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            production = directory / "eda-production.mdp"
            first = directory / "eda-nvt.mdp"
            second = directory / "thf-nvt.mdp"
            with patch("solvelec.cli.repository_root", return_value=ROOT):
                self.assertEqual(
                    main(
                        [
                            "render-gromacs-mdp",
                            "--protocol",
                            "pilot",
                            "--stage",
                            "production",
                            "--system",
                            "eda_1p5m",
                            "--replica",
                            "1",
                            "--output",
                            str(production),
                        ]
                    ),
                    0,
                )
                for system, output in (("eda_1p5m", first), ("pure_thf", second)):
                    result = main(
                        [
                            "render-gromacs-mdp",
                            "--protocol",
                            "pilot",
                            "--stage",
                            "nvt",
                            "--system",
                            system,
                            "--replica",
                            "1",
                            "--output",
                            str(output),
                        ]
                    )
                    self.assertEqual(result, 0)

            production_text = production.read_text(encoding="utf-8")
            eda_text = first.read_text(encoding="utf-8")
            thf_text = second.read_text(encoding="utf-8")
            self.assertIn("nsteps                      = 10000000", production_text)
            self.assertIn("dt                          = 0.002000", production_text)
            self.assertIn("nstxout-compressed          = 5000", production_text)
            self.assertNotEqual(eda_text, thf_text)

    def test_written_spec_records_the_initial_box(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "spec.json"
            with patch("solvelec.cli.repository_root", return_value=ROOT):
                result = main(
                    [
                        "write-spec",
                        "--system",
                        "pure_thf",
                        "--replica",
                        "1",
                        "--output",
                        str(output),
                    ]
                )
            self.assertEqual(result, 0)
            specification = json.loads(output.read_text(encoding="utf-8"))
            self.assertGreater(specification["initial_box_angstrom"], 20.0)


if __name__ == "__main__":
    unittest.main()
