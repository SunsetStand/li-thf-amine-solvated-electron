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

    def test_production_gromacs_command_resumes_from_internal_checkpoint(self) -> None:
        module = load_script("run_gromacs_stage")
        _grompp, mdrun = module.build_commands(
            "production",
            Path("production.mdp"),
            Path("npt.gro"),
            Path("topol.top"),
            Path("production"),
            4,
            Path("npt.cpt"),
            Path(".resume/production.cpt"),
            10,
        )
        self.assertEqual(Path(mdrun[mdrun.index("-cpi") + 1]), Path(".resume/production.cpt"))
        self.assertIn("-append", mdrun)
        self.assertEqual(mdrun[mdrun.index("-cpt") + 1], "10")

    def test_restart_workspace_is_reused_only_for_matching_inputs(self) -> None:
        module = load_script("run_gromacs_stage")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = []
            for name in ("nvt.mdp", "em.gro", "topol.top"):
                path = root / name
                path.write_text(name, encoding="utf-8")
                inputs.append(path)
            fingerprint = module.input_fingerprint(
                "nvt", inputs[0], inputs[1], inputs[2], None, 4, 15
            )
            workspace, resumed = module.prepare_restart_workspace(root / "output", fingerprint)
            self.assertFalse(resumed)
            (workspace / "nvt.tpr").write_text("tpr", encoding="utf-8")
            (workspace / "nvt.cpt").write_text("cpt", encoding="utf-8")

            same_workspace, resumed = module.prepare_restart_workspace(root / "output", fingerprint)
            self.assertTrue(resumed)
            self.assertEqual(same_workspace, workspace)

            inputs[0].write_text("changed", encoding="utf-8")
            changed = module.input_fingerprint("nvt", inputs[0], inputs[1], inputs[2], None, 4, 15)
            _, resumed = module.prepare_restart_workspace(root / "output", changed)
            self.assertFalse(resumed)
            self.assertFalse((workspace / "nvt.cpt").exists())

    def test_restart_outputs_are_promoted_only_with_manifest_and_logs(self) -> None:
        module = load_script("run_gromacs_stage")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "output"
            workspace = output / ".resume"
            workspace.mkdir(parents=True)
            for name in (
                "nvt.tpr",
                "nvt.gro",
                "nvt.edr",
                "nvt.log",
                "nvt.cpt",
                "nvt.xtc",
                "grompp.log",
                "mdrun.stdout.log",
                "manifest.json",
            ):
                (workspace / name).write_text(name, encoding="utf-8")

            module.promote_outputs(
                workspace, output, "nvt", ["tpr", "gro", "edr", "log", "cpt", "xtc"]
            )

            self.assertFalse(workspace.exists())
            self.assertTrue((output / "nvt.xtc").is_file())
            self.assertTrue((output / "manifest.json").is_file())

    def test_classical_pilot_metrics_and_replica_gate(self) -> None:
        module = load_script("validate_classical_pilot")
        metrics = module.replica_metrics(
            total_mass_g_mol=64 * 72.107,
            volumes_nm3=[8.80, 8.82, 8.81, 8.79, 8.80],
            times_ps=[0.0, 5000.0, 10000.0, 15000.0, 20000.0],
            amine_count=0,
            target_concentration_m=0.0,
            expected_duration_ns=20.0,
            concentration_tolerance_m=0.05,
            minimum_trajectory_fraction=0.98,
            density_half_relative_tolerance=0.02,
            engine_converged=True,
        )
        self.assertTrue(metrics["ready"])
        self.assertAlmostEqual(metrics["sampled_duration_ns"], 20.0)
        self.assertGreater(metrics["mean_density_g_ml"], 0.8)

        records = [
            {
                "system_id": "pure_thf",
                "replica": replica,
                "metrics": {**metrics, "mean_density_g_ml": density},
            }
            for replica, density in enumerate((0.87, 0.875, 0.872), start=1)
        ]
        summary = module.summarize_records(records, 0.03)
        self.assertTrue(summary["ready"])
        self.assertEqual(summary["systems"]["pure_thf"]["replicas"], [1, 2, 3])

        records[0]["metrics"] = {"ready": False, "error": "trajectory unreadable"}
        failed_summary = module.summarize_records(records, 0.03)
        self.assertFalse(failed_summary["ready"])
        self.assertFalse(
            failed_summary["systems"]["pure_thf"]["checks"]["replica_density_consistent"]
        )


if __name__ == "__main__":
    unittest.main()
