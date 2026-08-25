from __future__ import annotations

import os
import stat
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from solvelec.cli import main

ROOT = Path(__file__).resolve().parents[1]
RUN_SH = ROOT / "run.sh"
SLURM_DRIVER = ROOT / "configs" / "slurm" / "tmc-amd-driver.sbatch"
STAGE_MODULES = ROOT / "configs" / "slurm" / "tmc-amd-stage-modules.sh"
STORAGE_HELPER = ROOT / "configs" / "slurm" / "tmc-amd-storage.sh"
TMC_PROFILE = ROOT / "configs" / "profiles" / "tmc-amd" / "config.v9+.yaml"
ENGINE_RUNNER = ROOT / "workflow" / "scripts" / "run_checked_engine.py"
ENVIRONMENT_DIR = ROOT / "configs" / "environments"


class SlurmSafetyTests(unittest.TestCase):
    def test_shell_entry_points_have_valid_bash_syntax(self) -> None:
        for script in (RUN_SH, SLURM_DRIVER, STAGE_MODULES, STORAGE_HELPER):
            with self.subTest(script=script):
                completed = subprocess.run(
                    ["bash", "-n", str(script)], capture_output=True, text=True, check=False
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_tmc_profile_enforces_partition_and_four_cpu_minimum(self) -> None:
        profile = TMC_PROFILE.read_text(encoding="utf-8")
        self.assertIn("slurm_partition=amd", profile)
        self.assertIn("cpus_per_task=4", profile)
        self.assertIn("slurm-no-account: true", profile)

    def test_tmc_engine_stages_pin_separate_mpi_module_families(self) -> None:
        modules = STAGE_MODULES.read_text(encoding="utf-8")
        self.assertIn("module load cp2k/2023.2", modules)
        self.assertIn("module load gromacs/2023", modules)
        self.assertIn("module load orca/6.1.1", modules)
        driver = SLURM_DRIVER.read_text(encoding="utf-8")
        self.assertIn("incompatible MPI modules", driver)

    def test_tmc_storage_is_scoped_below_the_current_user(self) -> None:
        helper = STORAGE_HELPER.read_text(encoding="utf-8")
        self.assertIn("/data/home/storage/Backup_Data/${user_name}/", helper)
        self.assertIn("Rejected storage root", helper)
        self.assertIn('"${root}/runs"', helper)
        self.assertIn('"${root}/software/conda/envs"', helper)

    def test_tmc_tool_environments_are_split_at_mpi_boundaries(self) -> None:
        ambertools = (ENVIRONMENT_DIR / "tmc-ambertools.yml").read_text(encoding="utf-8")
        chem = (ENVIRONMENT_DIR / "tmc-chem-tools.yml").read_text(encoding="utf-8")
        qe = (ENVIRONMENT_DIR / "tmc-qe.yml").read_text(encoding="utf-8")
        self.assertIn("ambertools=26.0=cuda_None_nompi_py312*", ambertools)
        self.assertNotIn("packmol", ambertools)
        for package in ("packmol", "openbabel", "xtb", "crest"):
            self.assertIn(package, chem)
        self.assertIn("qe>=7,<8", qe)
        self.assertIn("openmpi>=4,<5", qe)

        stages = STAGE_MODULES.read_text(encoding="utf-8")
        self.assertLess(
            stages.index("solvelec_activate_support_tools"),
            stages.index("module load gromacs/2023"),
        )
        self.assertIn("solvelec_activate_qe", stages)

    def test_storage_and_tool_install_commands_require_slurm(self) -> None:
        runner = RUN_SH.read_text(encoding="utf-8")
        for command in ("storage-init", "tools-install"):
            self.assertIn(f"{command} must run inside a Slurm allocation", runner)
        self.assertIn("module load miniconda3", runner)
        self.assertNotIn(
            'if [[ "$1" == "tools-install" ]]',
            SLURM_DRIVER.read_text(encoding="utf-8"),
        )

    def test_logs_command_is_lightweight_and_available_without_python(self) -> None:
        environment = os.environ.copy()
        environment.pop("SLURM_JOB_ID", None)
        environment.pop("SOLVELEC_REQUIRE_SLURM", None)
        completed = subprocess.run(
            ["bash", str(RUN_SH), "logs", "probe"],
            cwd=ROOT,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    @unittest.skipIf(os.name == "nt", "fake executable routing is covered by Linux CI")
    def test_task_command_is_submitted_when_sbatch_is_available(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            fake_bin = temporary / "bin"
            fake_bin.mkdir()
            capture = temporary / "sbatch-arguments.txt"
            sbatch = fake_bin / "sbatch"
            sbatch.write_text(
                "#!/usr/bin/env bash\n"
                'printf \'%s\\n\' "$@" > "$SOLVELEC_CAPTURE"\n'
                "printf '12345\\n'\n",
                encoding="utf-8",
            )
            sbatch.chmod(sbatch.stat().st_mode | stat.S_IXUSR)
            environment = os.environ.copy()
            environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
            environment["SOLVELEC_CAPTURE"] = str(capture)
            environment.pop("SLURM_JOB_ID", None)

            completed = subprocess.run(
                ["bash", str(RUN_SH), "doctor", "--require", "input_bundle"],
                cwd=ROOT,
                env=environment,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("Submitted doctor as Slurm job 12345", completed.stdout)
            arguments = capture.read_text(encoding="utf-8").splitlines()
            self.assertIn("--partition=amd", arguments)
            self.assertIn("--cpus-per-task=4", arguments)
            self.assertIn(f"--export=ALL,SOLVELEC_REQUIRE_SLURM=1,SOLVELEC_ROOT={ROOT}", arguments)
            self.assertIn("doctor", arguments)
            self.assertIn("input_bundle", arguments)

    def test_driver_uses_exported_or_slurm_submission_directory(self) -> None:
        driver = SLURM_DRIVER.read_text(encoding="utf-8")
        self.assertIn("SOLVELEC_ROOT", driver)
        self.assertIn("SLURM_SUBMIT_DIR", driver)
        self.assertNotIn('dirname "${BASH_SOURCE[0]}"', driver)

    def test_engine_runner_refuses_outside_required_slurm_allocation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "engine.out"
            environment = os.environ.copy()
            environment["SOLVELEC_REQUIRE_SLURM"] = "1"
            environment.pop("SLURM_JOB_ID", None)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(ENGINE_RUNNER),
                    "--engine",
                    "cp2k",
                    "--output",
                    str(output),
                    "--",
                    "command-that-must-not-run",
                ],
                cwd=ROOT,
                env=environment,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 2)
            self.assertIn("outside a Slurm allocation", completed.stderr)
            self.assertFalse(output.exists())

    def test_direct_python_cli_refuses_when_sbatch_is_visible(self) -> None:
        environment = {"SOLVELEC_REQUIRE_SLURM": "0"}
        with (
            patch.dict(os.environ, environment, clear=True),
            patch("solvelec.cli.shutil.which", return_value="/usr/bin/sbatch"),
            redirect_stderr(StringIO()) as stderr,
        ):
            return_code = main(["validate"])

        self.assertEqual(return_code, 2)
        self.assertIn("use ./run.sh instead", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
