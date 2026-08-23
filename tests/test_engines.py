from __future__ import annotations

import subprocess
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from solvelec.cli import main
from solvelec.engines import (
    EngineSpec,
    EngineStatus,
    _find_executable,
    detect_engine,
    doctor_report,
    required_engine_names,
)


class EngineDetectionTests(unittest.TestCase):
    def test_finds_executable_beside_active_venv_python(self) -> None:
        spec = EngineSpec("snakemake", ("snakemake",), category="workflow")
        python_executable = "/repo/.venv/bin/python"
        python_bin = str(Path(python_executable).parent)

        def fake_which(name: str, path: str | None = None) -> str | None:
            if name == "snakemake" and path == python_bin:
                return "/repo/.venv/bin/snakemake"
            return None

        with (
            patch("solvelec.engines.sys.executable", python_executable),
            patch("solvelec.engines.shutil.which", side_effect=fake_which),
        ):
            self.assertEqual(_find_executable(spec), "/repo/.venv/bin/snakemake")

    def test_rejects_linux_orca_screen_reader(self) -> None:
        spec = EngineSpec(
            "orca",
            ("orca",),
            version_args=(),
            accepted_output_markers=("program version", "o   r   c   a"),
            rejected_output_markers=("screen reader",),
        )
        completed = subprocess.CompletedProcess(
            ["/usr/bin/orca"],
            1,
            stdout="",
            stderr="Cannot start the screen reader because it cannot connect to the Desktop.\n",
        )
        with (
            patch("solvelec.engines.shutil.which", return_value="/usr/bin/orca"),
            patch("solvelec.engines.subprocess.run", return_value=completed),
        ):
            status = detect_engine(spec)

        self.assertFalse(status.found)
        self.assertEqual(status.executable, "/usr/bin/orca")
        self.assertIn("rejected marker", status.error or "")

    def test_accepts_quantum_chemistry_orca_banner(self) -> None:
        spec = EngineSpec(
            "orca",
            ("orca",),
            version_args=(),
            accepted_output_markers=("program version", "o   r   c   a"),
            rejected_output_markers=("screen reader",),
        )
        completed = subprocess.CompletedProcess(
            ["/opt/orca/orca"],
            0,
            stdout="O   R   C   A\nProgram Version 6.0.1\n",
            stderr="",
        )
        with (
            patch("solvelec.engines.shutil.which", return_value="/opt/orca/orca"),
            patch("solvelec.engines.subprocess.run", return_value=completed),
        ):
            status = detect_engine(spec)

        self.assertTrue(status.found)
        self.assertIsNone(status.error)

    def test_accepts_orca_611_missing_parameterfile_signature(self) -> None:
        spec = EngineSpec(
            "orca",
            ("orca",),
            version_args=(),
            accepted_output_markers=(
                "program version",
                "o   r   c   a",
                "requires the name of a parameterfile as argument",
            ),
            rejected_output_markers=("screen reader",),
        )
        completed = subprocess.CompletedProcess(
            ["/opt/orca/orca"],
            255,
            stdout="This program requires the name of a parameterfile as argument\n",
            stderr="",
        )
        with (
            patch("solvelec.engines.shutil.which", return_value="/opt/orca/orca"),
            patch("solvelec.engines.subprocess.run", return_value=completed),
        ):
            status = detect_engine(spec)

        self.assertTrue(status.found)
        self.assertIsNone(status.error)

    def test_engine_specific_timeout_is_used(self) -> None:
        spec = EngineSpec("slow", ("slow",), timeout_seconds=30.0)
        completed = subprocess.CompletedProcess(["/opt/slow", "--version"], 0, "slow 1\n", "")
        with (
            patch("solvelec.engines.shutil.which", return_value="/opt/slow"),
            patch("solvelec.engines.subprocess.run", return_value=completed) as run,
        ):
            status = detect_engine(spec)

        self.assertTrue(status.found)
        self.assertEqual(run.call_args.kwargs["timeout"], 30.0)


class DoctorRequirementTests(unittest.TestCase):
    def test_input_bundle_requires_only_git_and_snakemake(self) -> None:
        self.assertEqual(required_engine_names(["input_bundle"]), frozenset({"git", "snakemake"}))

    def test_requirements_can_be_combined(self) -> None:
        required = required_engine_names(["cdft", "hpc"])
        self.assertEqual(required, frozenset({"git", "snakemake", "cp2k", "slurm", "mpi"}))

    def test_mpi_is_part_of_parallel_engine_stage_gates(self) -> None:
        for stage in ("classical_md", "cdft", "embedded_vde", "plane_wave_gate"):
            with self.subTest(stage=stage):
                self.assertIn("mpi", required_engine_names([stage]))

    def test_unknown_requirement_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "unknown doctor requirement"):
            required_engine_names(["not-a-stage"])

    def test_report_distinguishes_inventory_from_required_failures(self) -> None:
        specs = (
            EngineSpec("git", ("git",), category="core"),
            EngineSpec("snakemake", ("snakemake",), category="workflow"),
            EngineSpec("orca", ("orca",)),
        )
        statuses = {
            "git": EngineStatus("git", True, "/usr/bin/git", "git 2", "core", "all"),
            "snakemake": EngineStatus(
                "snakemake", False, None, None, "workflow", "workflow", "not found"
            ),
            "orca": EngineStatus(
                "orca", False, "/usr/bin/orca", "screen reader", "chemistry", "embedded_vde"
            ),
        }
        with patch("solvelec.engines.detect_engine", side_effect=lambda spec: statuses[spec.name]):
            report = doctor_report(specs, requirements=["input_bundle"])

        self.assertFalse(report["ready"])
        self.assertEqual(report["missing_required"], ["snakemake"])
        self.assertIn("orca", report["missing"])
        engine_records = {item["name"]: item for item in report["engines"]}
        self.assertTrue(engine_records["snakemake"]["required"])
        self.assertFalse(engine_records["orca"]["required"])

    def test_cli_strict_engines_remains_production_alias(self) -> None:
        report = {"ready": False, "engines": [], "missing": [], "missing_required": ["orca"]}
        with (
            patch("solvelec.cli.doctor_report", return_value=report) as mocked,
            redirect_stdout(StringIO()),
        ):
            return_code = main(["doctor", "--strict-engines"])

        self.assertEqual(return_code, 3)
        mocked.assert_called_once_with(requirements=["production"])


if __name__ == "__main__":
    unittest.main()
