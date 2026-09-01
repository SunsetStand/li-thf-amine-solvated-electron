from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import ModuleType
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str) -> ModuleType:
    path = ROOT / "workflow" / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ClassicalWorkflowTests(unittest.TestCase):
    def test_mol2_summary_counts_atoms_and_charge(self) -> None:
        module = load_script("prepare_molecule")
        with tempfile.TemporaryDirectory() as directory:
            mol2 = Path(directory) / "molecule.mol2"
            mol2.write_text(
                "@<TRIPOS>MOLECULE\nTHF\n"
                "@<TRIPOS>ATOM\n"
                "1 C1 0 0 0 c3 1 THF 0.125\n"
                "2 O1 0 0 0 os 1 THF -0.125 STATUS\n"
                "@<TRIPOS>BOND\n",
                encoding="utf-8",
            )
            self.assertEqual(module.mol2_summary(mol2), (2, 0.0))

    def test_parameterization_command_uses_isolated_working_directory(self) -> None:
        module = load_script("prepare_molecule")
        with tempfile.TemporaryDirectory() as directory:
            working_directory = Path(directory) / "scratch"
            working_directory.mkdir()
            log = Path(directory) / "parameterize.log"
            with patch.object(module.subprocess, "run") as run:
                run.return_value.returncode = 0
                module.run_checked(["antechamber", "-h"], log, working_directory)

            self.assertEqual(run.call_args.kwargs["cwd"], working_directory)

    def test_tleap_input_registers_gaff2_templates_and_box(self) -> None:
        module = load_script("build_gromacs_system")
        text = module.render_leap_input(
            Path("packed.pdb"),
            22.5,
            [("THF", Path("thf.mol2"), Path("thf.frcmod"), 64)],
            Path("system.prmtop"),
            Path("system.inpcrd"),
        )
        self.assertIn("source leaprc.gaff2", text)
        self.assertIn('loadamberparams "thf.frcmod"', text)
        self.assertIn('THF = loadmol2 "thf.mol2"', text)
        self.assertIn("check THF", text)
        self.assertNotIn("MOL0", text)
        self.assertIn("set system box { 22.50000000 22.50000000 22.50000000 }", text)
        self.assertIn('saveamberparm system "system.prmtop" "system.inpcrd"', text)

    def test_tleap_input_rejects_invalid_or_duplicate_template_names(self) -> None:
        module = load_script("build_gromacs_system")
        common = (Path("packed.pdb"), 22.5)
        outputs = (Path("system.prmtop"), Path("system.inpcrd"))

        with self.assertRaisesRegex(ValueError, "valid TLeap identifier"):
            module.render_leap_input(
                *common,
                [("1BAD", Path("bad.mol2"), Path("bad.frcmod"), 1)],
                *outputs,
            )
        with self.assertRaisesRegex(ValueError, "duplicate TLeap residue template"):
            module.render_leap_input(
                *common,
                [
                    ("THF", Path("thf.mol2"), Path("thf.frcmod"), 32),
                    ("THF", Path("thf-2.mol2"), Path("thf-2.frcmod"), 32),
                ],
                *outputs,
            )

    def test_npt_gromacs_command_uses_checkpoint_and_four_threads(self) -> None:
        module = load_script("run_gromacs_stage")
        grompp, mdrun = module.build_commands(
            "npt",
            Path("npt.mdp"),
            Path("nvt.gro"),
            Path("topol.top"),
            Path("npt"),
            4,
            Path("nvt.cpt"),
        )
        self.assertEqual(grompp[grompp.index("-t") + 1], "nvt.cpt")
        self.assertEqual(mdrun[mdrun.index("-ntomp") + 1], "4")
        self.assertEqual(mdrun[mdrun.index("-ntmpi") + 1], "1")


if __name__ == "__main__":
    unittest.main()
