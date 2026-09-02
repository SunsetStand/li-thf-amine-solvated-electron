from __future__ import annotations

import os
import shutil
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
SNAKEFILE = ROOT / "workflow" / "Snakefile"
ENGINE_RUNNER = ROOT / "workflow" / "scripts" / "run_checked_engine.py"
STAGE_RUNNER = ROOT / "workflow" / "scripts" / "run_tmc_stage.sh"
ENVIRONMENT_DIR = ROOT / "configs" / "environments"


class SlurmSafetyTests(unittest.TestCase):
    def test_shell_entry_points_have_valid_bash_syntax(self) -> None:
        for script in (RUN_SH, SLURM_DRIVER, STAGE_MODULES, STORAGE_HELPER, STAGE_RUNNER):
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
        self.assertIn("slurm-status-command: squeue", profile)

    def test_classical_pilot_is_an_explicit_restricted_target(self) -> None:
        runner = RUN_SH.read_text(encoding="utf-8")
        snakefile = SNAKEFILE.read_text(encoding="utf-8")
        self.assertIn("input_bundle|classical_smoke|classical_pilot", runner)
        self.assertIn("classical_analysis|snapshot_bank", runner)
        self.assertIn('if CAMPAIGN != "pilot"', snakefile)
        self.assertIn("gromacs_pilot_md:", TMC_PROFILE.read_text(encoding="utf-8"))
        self.assertIn("analyze_classical_replica:", TMC_PROFILE.read_text(encoding="utf-8"))

    def test_stage_a_cannot_reschedule_completed_classical_md(self) -> None:
        rules = (ROOT / "workflow" / "rules" / "40_classical_analysis.smk").read_text(
            encoding="utf-8"
        )
        analyze_rule, snapshot_rules = rules.split("rule summarize_classical_analysis:", 1)
        self.assertIn("tpr=lambda wildcards", analyze_rule)
        self.assertIn("trajectory=lambda wildcards", analyze_rule)
        self.assertNotIn("tpr=(", analyze_rule)
        self.assertNotIn("trajectory=(", analyze_rule)
        snapshot_rule = snapshot_rules.split("rule select_classical_snapshot:", 1)[1]
        self.assertIn("tpr=lambda wildcards", snapshot_rule)
        self.assertIn("trajectory=lambda wildcards", snapshot_rule)

    def test_stage_b_treats_stage_a_as_an_immutable_handoff(self) -> None:
        rules = (ROOT / "workflow" / "rules" / "55_stage_b.smk").read_text(
            encoding="utf-8"
        )
        prepare_rule = rules.split("rule summarize_stage_b_candidates:", 1)[0]
        for name in ("spec", "snapshot_metadata", "xyz", "cell"):
            self.assertIn(f"{name}=lambda wildcards", prepare_rule)
        self.assertNotIn("snapshot_bank.done", rules)
        runner = RUN_SH.read_text(encoding="utf-8")
        self.assertIn("stage_b_candidates|stage_b", runner)
        profile = TMC_PROFILE.read_text(encoding="utf-8")
        self.assertIn("run_stage_b_cp2k_smoke:", profile)
        self.assertIn("tasks: 32", profile)
        self.assertIn("cpus_per_task: 1", profile)

    def test_slurm_plugin_has_compute_node_hang_fix(self) -> None:
        project = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        pixi = (ROOT / "pixi.toml").read_text(encoding="utf-8")
        self.assertIn("snakemake-executor-plugin-slurm>=2.7,<3", project)
        self.assertIn('snakemake-executor-plugin-slurm = ">=2.7,<3"', pixi)

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
        self.assertIn("CONDA_ALWAYS_YES=true", runner)
        self.assertNotIn("conda env update --yes", runner)
        self.assertNotIn("conda env create --yes", runner)
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

    def test_nested_slurm_controller_does_not_inherit_parent_job_context(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            temporary_run = temporary / "run.sh"
            shutil.copy2(RUN_SH, temporary_run)

            slurm_config = temporary / "configs" / "slurm"
            profile = temporary / "configs" / "profiles" / "tmc-amd"
            fake_bin = temporary / ".venv" / "bin"
            slurm_config.mkdir(parents=True)
            profile.mkdir(parents=True)
            fake_bin.mkdir(parents=True)
            shutil.copy2(STORAGE_HELPER, slurm_config / STORAGE_HELPER.name)
            (profile / "config.v9+.yaml").write_text("executor: slurm\n", encoding="utf-8")

            python = fake_bin / "python"
            python.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
            python.chmod(python.stat().st_mode | stat.S_IXUSR)

            capture = temporary / "snakemake-environment.txt"
            snakemake = fake_bin / "snakemake"
            snakemake.write_text(
                '#!/usr/bin/env bash\nenv | sort > "$SOLVELEC_CAPTURE"\n',
                encoding="utf-8",
            )
            snakemake.chmod(snakemake.stat().st_mode | stat.S_IXUSR)

            environment = os.environ.copy()
            environment["SLURM_JOB_ID"] = "parent-123"
            environment["SLURM_JOB_NAME"] = "solvelec-submit"
            environment["SLURM_CPUS_PER_TASK"] = "4"
            environment["SOLVELEC_DEFAULT_PROFILE"] = "tmc-amd"
            environment["SOLVELEC_REQUIRE_SLURM"] = "1"
            environment["SOLVELEC_CAPTURE"] = str(capture)
            environment["USER"] = "solvelec-test"
            environment["SOLVELEC_STORAGE_ROOT"] = (
                "/data/home/storage/Backup_Data/solvelec-test/project"
            )

            completed = subprocess.run(
                [
                    "bash",
                    str(temporary_run),
                    "submit",
                    "--campaign",
                    "smoke",
                    "--target",
                    "classical_smoke",
                ],
                cwd=temporary,
                env=environment,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            controller_environment = capture.read_text(encoding="utf-8").splitlines()
            self.assertIn("Workflow rules remain protected", completed.stdout)
            self.assertFalse(
                any(item.startswith("SLURM_") for item in controller_environment),
                f"{completed.stdout}\n{controller_environment}",
            )
            self.assertIn("SOLVELEC_REQUIRE_SLURM=1", controller_environment)

    def test_driver_uses_exported_or_slurm_submission_directory(self) -> None:
        driver = SLURM_DRIVER.read_text(encoding="utf-8")
        self.assertIn("SOLVELEC_ROOT", driver)
        self.assertIn("SLURM_SUBMIT_DIR", driver)
        self.assertNotIn('dirname "${BASH_SOURCE[0]}"', driver)

    def test_unlock_checks_for_active_controllers_before_clearing_lock(self) -> None:
        runner = RUN_SH.read_text(encoding="utf-8")
        queue_check = 'squeue -h -u "${USER}" -n solvelec-submit,solvelec-resume'
        self.assertIn(queue_check, runner)
        self.assertLess(
            runner.index(queue_check), runner.index('storage_root="${storage_root}" --unlock')
        )
        self.assertIn("refusing to unlock while workflow controllers are active", runner)

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

    def test_stage_runner_refuses_outside_slurm(self) -> None:
        environment = os.environ.copy()
        environment.pop("SLURM_JOB_ID", None)
        completed = subprocess.run(
            ["bash", str(STAGE_RUNNER), "classical_md", "--", "true"],
            cwd=ROOT,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("outside a Slurm allocation", completed.stderr)

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
