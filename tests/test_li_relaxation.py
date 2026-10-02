from __future__ import annotations

import copy
import csv
import importlib.util
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from solvelec.candidates import read_xyz, write_xyz
from solvelec.config import load_config
from solvelec.li_relaxation import (
    PAIRS,
    energy_metrics,
    last_xyz_frame,
    read_json,
    render_input,
    seed_li,
    settings_from,
    validate_output,
    write_json,
)
from solvelec.parsers import HARTREE_TO_EV

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "workflow/templates/cp2k/li_relaxation.inp.tpl"
METHODS = load_config(ROOT / "configs/methods.yaml")


def load_script():
    spec = importlib.util.spec_from_file_location(
        "li_relaxation_workflow", ROOT / "workflow/scripts/li_relaxation.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def output_text(job, energy=-10.0, *, deviation=0.0002, strength=0.3, converged=True):
    target = job["target_electrons"]
    total = energy + deviation * strength
    electrons = job["expected_valence_electrons"]
    geometry = "GEOMETRY OPTIMIZATION COMPLETED" if converged else ""
    return f""" Number of electrons: {(electrons + 1) // 2}
 Number of electrons: {(electrons - 1) // 2}
 CDFT SCF iter = 4 RMS gradient = {abs(deviation):.12e} energy = {total:.12f}
 CDFT SCF loop converged in 4 iterations or 20 steps
 Target value of constraint : {target:.12f}
 Current value of constraint : {target + deviation:.12f}
 Deviation from target : {deviation:.12e}
 Strength of constraint : {strength:.12f}
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: {total:.12f}
 {geometry}
 PROGRAM ENDED AT 2026-10-02
"""


class LiRelaxationTests(unittest.TestCase):
    def setUp(self):
        self.module = load_script()
        self.settings = settings_from(copy.deepcopy(METHODS))

    def make_bank(self, root, system="eda_3m"):
        methods = root / "methods.json"
        write_json(methods, METHODS)
        xyz = root / "seed.xyz"
        write_xyz(xyz, ["LI", "H", "H"], np.array([[2, 2, 2], [8, 8, 8], [8, 8, 9]]), "seed")
        cell = root / "cell.inc"
        cell.write_text("&CELL\n A 25 0 0\n B 0 25 0\n C 0 0 25\n PERIODIC XYZ\n&END CELL\n")
        bank = root / "bank.json"
        write_json(
            bank,
            {
                "ready": True,
                "system_id": system,
                "seeds": [
                    {
                        "snapshot": "s01",
                        "structure": {
                            "xyz": self.module.record(xyz),
                            "cell": self.module.record(cell),
                        },
                    }
                ],
            },
        )
        return methods, bank, xyz

    def render_job(
        self, root, methods, bank, coordinates, state, geometry, optimize, validation=None
    ):
        output = root / "cp2k.inp"
        self.module.render(
            SimpleNamespace(
                methods=str(methods),
                bank=str(bank),
                snapshot="s01",
                coordinates=str(coordinates),
                state=state,
                geometry_state=geometry,
                optimize=optimize,
                variant="pbe0_admm",
                template=str(TEMPLATE),
                output=str(output),
                geometry_validation=validation,
            )
        )
        return output.with_suffix(".job.json")

    def mock_data(self, root):
        data = root / "data"
        data.mkdir()
        (data / "BASIS_MOLOPT").write_text("H TZV2P-MOLOPT-GTH\nLi DZVP-MOLOPT-SR-GTH\n")
        (data / "BASIS_ADMM").write_text("H cFIT3\nLi cFIT3\n")
        for name in ("GTH_POTENTIALS", "dftd3.dat", "t_c_g.dat"):
            (data / name).touch()
        return data

    def mock_engine(self, job, energy, *, failure=False):
        def execute(command, *, cwd, env, stdout, stderr):
            self.assertEqual(env["OMP_NUM_THREADS"], "1")
            self.assertEqual(command[:1], ["mock-cp2k"])
            self.assertEqual(command[-2], "-i")
            text = output_text(job, energy=energy, converged=not failure)
            if failure:
                text += "SCF run NOT converged\n"
            stdout.write(text)
            if job["run_type"] == "GEO_OPT":
                elements, coordinates, _ = read_xyz(job["coordinates"]["path"])
                if job["state"] == "b":
                    coordinates[0, 0] += 0.1
                write_xyz(
                    Path(cwd) / "positions.xyz",
                    elements,
                    coordinates,
                    f"i = 3, E = {energy + 0.00006:.12f}",
                )
            else:
                for label in ("ELECTRON_DENSITY", "SPIN_DENSITY"):
                    (Path(cwd) / f"li_relax-{label}-1_0.cube").write_text("mock density fixture\n")
            return SimpleNamespace(returncode=0)

        return execute

    def test_energy_identity_and_sign_use_matched_differences(self):
        metrics = energy_metrics({"aa": -10, "ba": -9.9, "bb": -9.95, "ab": -9.98})
        self.assertAlmostEqual(metrics["vertical_transfer_ev"], 0.1 * HARTREE_TO_EV)
        self.assertAlmostEqual(metrics["relaxed_transfer_ev"], 0.05 * HARTREE_TO_EV)
        self.assertAlmostEqual(
            metrics["vertical_transfer_ev"] - metrics["product_relaxation_stabilization_ev"],
            metrics["relaxed_transfer_ev"],
        )
        for energies in ({"aa": -10}, {"aa": float("nan"), "ba": -9, "bb": -9, "ab": -10}):
            with self.assertRaises(ValueError):
                energy_metrics(energies)

    def test_existing_md_handoff_prepares_three_seeds_and_preserves_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            methods = root / "methods.json"
            write_json(methods, METHODS)
            paths = self.module.source_paths(root, self.settings, "eda_3m")
            for path in paths.values():
                path.parent.mkdir(parents=True, exist_ok=True)
            write_json(paths["spec"], {"amine": "eda", "component_counts": {"eda": 1}})
            paths["tpr"].write_bytes(b"immutable tpr")
            paths["trajectory"].write_bytes(b"immutable xtc")
            write_json(
                paths["analysis"],
                {
                    "ready": True,
                    "system_id": "eda_3m",
                    "replica": 1,
                    "inputs": {
                        key: self.module.record(paths[key]) for key in ("tpr", "trajectory", "spec")
                    },
                    "metrics": {
                        "autocorrelation": {
                            "density_g_ml": {"integrated_autocorrelation_time_ps": 100},
                            "volume_nm3": {"integrated_autocorrelation_time_ps": 100},
                        }
                    },
                },
            )
            times = [11000, 15000, 19000]
            paths["timeseries"].write_text(
                "frame_index,time_ps,elapsed_ps,density_g_ml,volume_nm3\n"
                + "".join(
                    f"{i},{t},{t},{0.9 + i / 100},{15 + i / 100}\n" for i, t in enumerate(times)
                )
            )
            source_hashes = {key: self.module.record(path)["sha256"] for key, path in paths.items()}
            positions = np.array([[0.5, 5, 5], [24.5, 5, 5], [12, 12, 12], [12, 12, 13]])

            class Atoms(list):
                def unwrap(self, **kwargs):
                    return positions.copy()

            atoms = Atoms(
                [
                    SimpleNamespace(index=i, name=element, element=element)
                    for i, element in enumerate(("N", "N", "H", "H"))
                ]
            )
            frames = [SimpleNamespace(time=t, dimensions=[25, 25, 25, 90, 90, 90]) for t in times]
            universe = SimpleNamespace(
                atoms=atoms, residues=[SimpleNamespace(atoms=atoms)], trajectory=frames
            )
            mda = SimpleNamespace(Universe=lambda *args: universe)
            args = SimpleNamespace(
                methods=str(methods),
                run_root=str(root),
                system="eda_3m",
                output_dir=str(root / "seeds"),
            )
            with patch.dict(sys.modules, {"MDAnalysis": mda}):
                self.assertEqual(self.module.prepare(args), 0)
            bank = read_json(root / "seeds/bank.json")
            self.assertEqual(bank["snapshot_count"], 3)
            self.assertEqual([seed["frame_index"] for seed in bank["seeds"]], [0, 1, 2])
            for seed in bank["seeds"]:
                elements, _, _ = read_xyz(seed["structure"]["xyz"]["path"])
                self.assertEqual(elements[0], "LI")
                self.assertGreater(seed["li_concentration_m"], 0)
            self.assertEqual(
                source_hashes,
                {key: self.module.record(path)["sha256"] for key, path in paths.items()},
            )
            paths["tpr"].write_bytes(b"different tpr")
            with (
                patch.dict(sys.modules, {"MDAnalysis": mda}),
                self.assertRaisesRegex(ValueError, "no longer matches"),
            ):
                self.module.prepare(args)

    def test_output_correction_and_convergence_electron_and_state_gates(self):
        job = {"target_electrons": 2, "expected_valence_electrons": 5, "run_type": "GEO_OPT"}
        text = output_text(job)
        value = validate_output(text, job, self.settings)
        self.assertTrue(value["ready"], value["problems"])
        self.assertAlmostEqual(value["physical_energy_hartree"], -10)
        self.assertAlmostEqual(value["constraint_energy_hartree"], 0.00006)
        failures = (
            text.replace("GEOMETRY OPTIMIZATION COMPLETED", ""),
            text + "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED\n",
            text + "CDFT SCF loop FAILED to converge\n",
            text.replace("Number of electrons: 3", "Number of electrons: 4"),
            output_text(job, deviation=0.02),
            text.replace(
                "Target value of constraint : 2.000", "Target value of constraint : 3.000"
            ),
            text.replace("energy [a.u.]: -9.999940000000", "energy [a.u.]: -9.888000000000"),
            text + "SCF run NOT converged\n",
        )
        for failed in failures:
            with self.subTest(output=failed):
                self.assertFalse(validate_output(failed, job, self.settings)["ready"])

    def test_hybrid_scaling_available_li_basis_and_cutoff_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, bank, xyz = self.make_bank(root)
            cell = read_json(bank)["seeds"][0]["structure"]["cell"]["path"]
            output = root / "input.inp"
            job = render_input(
                TEMPLATE,
                output,
                coordinates=xyz,
                cell=cell,
                settings=self.settings,
                state="b",
                variant="pbe0_admm",
                optimize=True,
            )
            text = output.read_text()
            self.assertIn("SCALE_X 0.75", text)
            self.assertIn("FRACTION 0.25", text)
            self.assertIn("BASIS_SET_FILE_NAME BASIS_ADMM\n", text)
            self.assertNotIn("BASIS_ADMM_MOLOPT", text)
            li_kind = text.split("&KIND Li")[1].split("&END KIND")[0]
            self.assertIn("BASIS_SET DZVP-MOLOPT-SR-GTH", li_kind)
            self.assertEqual((job["charge"], job["multiplicity"]), (0, 2))
            self.assertNotIn("GHOST", text)
            self.assertIn("CAVITY_CONFINE FALSE", text)
            Path(cell).write_text("A 9 0 0\nB 0 9 0\nC 0 0 9\n")
            with self.assertRaisesRegex(ValueError, "HFX cutoff"):
                render_input(
                    TEMPLATE,
                    output,
                    coordinates=xyz,
                    cell=cell,
                    settings=self.settings,
                    state="a",
                    variant="pbe0_admm",
                    optimize=False,
                )

    def test_periodic_li_seed_has_two_donors_and_clearance(self):
        elements = ["N", "N", "H"]
        positions = np.array([[0.5, 5, 5], [19.5, 5, 5], [10, 10, 10]])
        cell = np.eye(3) * 20
        li, record = seed_li(
            elements, positions, cell, [[0, 1]], distance_angstrom=2.2, minimum_scale=1
        )
        self.assertEqual(record["initial_coordination"], "bidentate")
        delta = positions[:2] - li
        delta -= 20 * np.rint(delta / 20)
        np.testing.assert_allclose(np.linalg.norm(delta, axis=1), [2.2, 2.2])
        self.assertGreaterEqual(record["minimum_contact_scale"], 1)
        with self.assertRaisesRegex(ValueError, "collision-free"):
            seed_li(elements, positions, cell, [[0, 1]], distance_angstrom=2.2, minimum_scale=100)

    def test_xyz_parser_rejects_truncated_final_frame(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "positions.xyz"
            path.write_text("1\nE = -1\nLi 0 0 0\n1\nE = -2\nLi 1 0 0\n")
            self.assertEqual(last_xyz_frame(path)[2], "E = -2")
            path.write_text(path.read_text() + "1\nE = -3\n")
            with self.assertRaisesRegex(ValueError, "truncated"):
                last_xyz_frame(path)

    def test_runner_rejects_no_slurm_and_preserves_failed_attempt(self):
        with patch.dict(os.environ, {}, clear=True), self.assertRaisesRegex(ValueError, "Slurm"):
            self.module.run(SimpleNamespace())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            methods, bank, xyz = self.make_bank(root)
            data = self.mock_data(root)
            work = root / "opt"
            job_path = self.render_job(work, methods, bank, xyz, "a", "seed", True)
            args = SimpleNamespace(
                job=str(job_path), workdir=str(work), engine_command=["--", "mock-cp2k"]
            )
            with patch.dict(os.environ, {"SLURM_JOB_ID": "123", "CP2K_DATA_DIR": str(data)}):
                with patch.object(
                    self.module.subprocess,
                    "run",
                    side_effect=self.mock_engine(read_json(job_path), -10, failure=True),
                ):
                    self.assertEqual(self.module.run(args), 4)
            self.assertTrue((work / "cp2k.out.failed").is_file())
            self.assertTrue((work / "validation.failed.json").is_file())
            self.assertFalse((work / "validation.json").exists())
            self.assertEqual(len(list(work.glob("attempts/*/cp2k.out"))), 1)

    def test_paired_optimization_to_four_point_analysis_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            methods, bank, xyz = self.make_bank(root)
            data = self.mock_data(root)
            # A local minimum need not be below a geometry in another basin.
            energies = {"aa": -10, "ba": -9.9, "bb": -9.95, "ab": -10.02}
            snapshot = root / "snapshot"
            geometry_validation = None
            with patch.dict(os.environ, {"SLURM_JOB_ID": "123", "CP2K_DATA_DIR": str(data)}):
                for state in ("a", "b"):
                    work = snapshot / "opt" / state
                    job_path = self.render_job(
                        work,
                        methods,
                        bank,
                        xyz,
                        state,
                        "seed" if state == "a" else "a",
                        True,
                        geometry_validation,
                    )
                    with patch.object(
                        self.module.subprocess,
                        "run",
                        side_effect=self.mock_engine(read_json(job_path), energies[state + state]),
                    ):
                        self.assertEqual(
                            self.module.run(
                                SimpleNamespace(
                                    job=str(job_path),
                                    workdir=str(work),
                                    engine_command=["--", "mock-cp2k"],
                                )
                            ),
                            0,
                        )
                    xyz, geometry_validation = work / "optimized.xyz", str(work / "validation.json")
                for pair in PAIRS:
                    geometry = snapshot / "opt" / pair[1]
                    work = snapshot / "sp/pbe0_admm" / pair
                    job_path = self.render_job(
                        work,
                        methods,
                        bank,
                        geometry / "optimized.xyz",
                        pair[0],
                        pair[1],
                        False,
                        str(geometry / "validation.json"),
                    )
                    with patch.object(
                        self.module.subprocess,
                        "run",
                        side_effect=self.mock_engine(read_json(job_path), energies[pair]),
                    ):
                        self.assertEqual(
                            self.module.run(
                                SimpleNamespace(
                                    job=str(job_path),
                                    workdir=str(work),
                                    engine_command=["--", "mock-cp2k"],
                                )
                            ),
                            0,
                        )
                    self.assertTrue((work / "spin.cube").exists())
                    self.assertTrue((work / "electron.cube").exists())
            output = root / "analysis.json"
            args = SimpleNamespace(
                methods=str(methods),
                bank=str(bank),
                snapshot="s01",
                variant="pbe0_admm",
                snapshot_dir=str(snapshot),
                output=str(output),
            )
            self.module.analyze(args)
            result = read_json(output)
            self.assertTrue(result["ready"], result["problems"])
            self.assertEqual(len(result["diagnostics"]), 1)
            self.assertAlmostEqual(
                result["metrics_ev"]["relaxed_transfer_ev"], 0.05 * HARTREE_TO_EV
            )
            validation_path = snapshot / "sp/pbe0_admm/ba/validation.json"
            original = read_json(validation_path)
            mixed = copy.deepcopy(original)
            mixed["job"]["geometry_state"] = "b"
            write_json(validation_path, mixed)
            with self.assertRaisesRegex(ValueError, "mixed electronic states"):
                self.module.analyze(args)
            mixed = copy.deepcopy(original)
            mixed["job"]["coordinates"] = self.module.record(snapshot / "opt/b/optimized.xyz")
            write_json(validation_path, mixed)
            with self.assertRaisesRegex(ValueError, "identical coordinates"):
                self.module.analyze(args)
            write_json(validation_path, original)
            changed = snapshot / "opt/a/optimized.xyz"
            changed.write_text(changed.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "checksum changed"):
                self.module.analyze(args)

    def test_summary_failed_record_keeps_partial_csv_and_gate_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            methods = root / "methods.json"
            write_json(methods, METHODS)
            paths = []
            for system, ready in (("eda_3m", True), ("tmeda_3m", False)):
                path = root / f"{system}.json"
                write_json(
                    path,
                    {
                        "system_id": system,
                        "snapshot": "s01",
                        "variant": "pbe0_admm",
                        "ready": ready,
                        "metrics_ev": {"relaxed_transfer_ev": 1} if ready else None,
                    },
                )
                paths.append(str(path))
            output = root / "summary.json"
            args = SimpleNamespace(
                methods=str(methods),
                records=paths,
                snapshots=["s01"],
                variants=["pbe0_admm"],
                campaign="li_eda_tmeda",
                output=str(output),
            )
            self.module.summary(args)
            summary = read_json(output)
            self.assertFalse(summary["ready"])
            self.assertEqual(summary["eda_minus_tmeda_ev"], {})
            with output.with_suffix(".csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["relaxed_transfer_ev"], "1")
            self.assertEqual(rows[1]["relaxed_transfer_ev"], "")
            with self.assertRaisesRegex(ValueError, "did not all pass"):
                self.module.gate(SimpleNamespace(summary=str(output), output=str(root / "done")))
            value = read_json(paths[1])
            value.update(ready=True, metrics_ev={"relaxed_transfer_ev": 2})
            write_json(paths[1], value)
            self.module.summary(args)
            summary = read_json(output)
            self.assertEqual(summary["eda_minus_tmeda_ev"]["pbe0_admm"]["relaxed_transfer_ev"], -1)
            stats = summary["systems_by_variant"]["pbe0_admm"]["eda_3m"]["metrics"]
            self.assertIsNone(stats["relaxed_transfer_ev"]["standard_deviation_ev"])
            self.module.gate(SimpleNamespace(summary=str(output), output=str(root / "done")))
            self.assertTrue((root / "done").exists())


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "optional workflow dependency")
class LiRelaxationDagTests(unittest.TestCase):
    def test_targets_reuse_source_files_without_md_rules(self):
        module = load_script()
        settings = settings_from(METHODS)
        with tempfile.TemporaryDirectory() as tmp:
            for system in settings["systems"]:
                for path in module.source_paths(tmp, settings, system).values():
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.touch()
            for target, optimizations, points in (
                ("li_relaxation_inputs", 0, 0),
                ("li_relaxation_pilot", 4, 8),
                ("li_relaxation", 12, 24),
                ("li_relaxation_benchmark", 4, 16),
            ):
                command = [
                    sys.executable,
                    "-m",
                    "snakemake",
                    "--snakefile",
                    "workflow/Snakefile",
                    "--directory",
                    str(ROOT),
                    target,
                    "--config",
                    "campaign=li_eda_tmeda",
                    f"run_root={tmp}",
                    "--dry-run",
                ]
                completed = subprocess.run(
                    command,
                    cwd=ROOT,
                    text=True,
                    capture_output=True,
                    env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
                    timeout=30,
                )
                self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
                table = completed.stdout.split("Job stats:", 1)[1].split("\n\n", 1)[0]
                rules = dict(re.findall(r"^([a-z_]+)\s+(\d+)\s*$", table, re.M))
                self.assertEqual(int(rules.get("run_li_relaxation_opt", 0)), optimizations)
                self.assertEqual(int(rules.get("run_li_relaxation_sp", 0)), points)
                self.assertTrue(all("li_relaxation" in rule or rule == "total" for rule in rules))
                if points:
                    self.assertIn("validate_li_relaxation_summary", rules)


if __name__ == "__main__":
    unittest.main()
