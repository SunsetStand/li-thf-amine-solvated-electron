from __future__ import annotations

import platform
import shutil
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class EngineSpec:
    name: str
    executables: tuple[str, ...]
    version_args: tuple[str, ...] = ("--version",)
    category: str = "chemistry"
    required_for: str = "production"


@dataclass(frozen=True)
class EngineStatus:
    name: str
    found: bool
    executable: str | None
    version: str | None
    category: str
    required_for: str
    error: str | None = None

    def as_dict(self) -> dict[str, str | bool | None]:
        return asdict(self)


ENGINE_SPECS: tuple[EngineSpec, ...] = (
    EngineSpec("git", ("git",), category="core", required_for="all"),
    EngineSpec("packmol", ("packmol",), version_args=("-h",), required_for="classical_md"),
    EngineSpec("openbabel", ("obabel",), required_for="molecule_generation"),
    EngineSpec("ambertools", ("antechamber",), version_args=("-h",), required_for="classical_md"),
    EngineSpec("gromacs", ("gmx", "gmx_mpi"), required_for="classical_md"),
    EngineSpec("cp2k", ("cp2k.psmp", "cp2k.popt", "cp2k"), required_for="cdft"),
    EngineSpec("orca", ("orca",), version_args=(), required_for="embedded_vde"),
    EngineSpec("quantum_espresso", ("pw.x",), version_args=("-h",), required_for="plane_wave_gate"),
    EngineSpec("xtb", ("xtb",), required_for="conformer_search"),
    EngineSpec("crest", ("crest",), required_for="conformer_search"),
    EngineSpec("snakemake", ("snakemake",), category="workflow", required_for="workflow"),
    EngineSpec(
        "apptainer", ("apptainer", "singularity"), category="workflow", required_for="container"
    ),
    EngineSpec(
        "slurm", ("sbatch",), version_args=("--version",), category="scheduler", required_for="hpc"
    ),
    EngineSpec(
        "mpi",
        ("mpirun", "srun"),
        version_args=("--version",),
        category="scheduler",
        required_for="hpc",
    ),
)


def _first_line(text: str) -> str | None:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped[:300]
    return None


def detect_engine(spec: EngineSpec, timeout_seconds: float = 5.0) -> EngineStatus:
    executable = next((shutil.which(name) for name in spec.executables if shutil.which(name)), None)
    if not executable:
        return EngineStatus(
            spec.name, False, None, None, spec.category, spec.required_for, "not found on PATH"
        )
    version = None
    error = None
    try:
        completed = subprocess.run(
            [executable, *spec.version_args],
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
        version = _first_line(completed.stdout) or _first_line(completed.stderr)
        if completed.returncode not in (0, 1) and version is None:
            error = f"version command exited {completed.returncode}"
    except (OSError, subprocess.SubprocessError) as exc:
        error = str(exc)
    return EngineStatus(
        spec.name, True, executable, version, spec.category, spec.required_for, error
    )


def doctor_report(specs: Sequence[EngineSpec] = ENGINE_SPECS) -> dict[str, object]:
    statuses = [detect_engine(spec) for spec in specs]
    return {
        "python": {"version": sys.version.split()[0], "executable": sys.executable},
        "platform": platform.platform(),
        "engines": [status.as_dict() for status in statuses],
        "missing": [status.name for status in statuses if not status.found],
    }
