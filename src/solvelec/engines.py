from __future__ import annotations

import platform
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class EngineSpec:
    name: str
    executables: tuple[str, ...]
    version_args: tuple[str, ...] = ("--version",)
    timeout_seconds: float = 15.0
    category: str = "chemistry"
    required_for: str = "production"
    accepted_output_markers: tuple[str, ...] = ()
    rejected_output_markers: tuple[str, ...] = ()


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
    EngineSpec(
        "cp2k",
        ("cp2k.psmp", "cp2k.popt", "cp2k"),
        timeout_seconds=30.0,
        required_for="cdft",
    ),
    EngineSpec(
        "orca",
        ("orca",),
        version_args=(),
        required_for="embedded_vde",
        accepted_output_markers=(
            "program version",
            "o   r   c   a",
            "requires the name of a parameterfile as argument",
        ),
        rejected_output_markers=("screen reader",),
    ),
    EngineSpec("quantum_espresso", ("pw.x",), version_args=("-h",), required_for="plane_wave_gate"),
    EngineSpec("xtb", ("xtb",), required_for="conformer_search"),
    EngineSpec("crest", ("crest",), required_for="conformer_search"),
    EngineSpec(
        "snakemake",
        ("snakemake",),
        timeout_seconds=30.0,
        category="workflow",
        required_for="workflow",
    ),
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


# These are workflow capability gates, not package-installation groups.  A user can
# therefore ask whether the next stage is runnable without being told that every
# optional production engine must already be installed.
REQUIREMENT_ENGINES: dict[str, tuple[str, ...]] = {
    "input_bundle": ("git", "snakemake"),
    "molecule_generation": ("git", "snakemake", "openbabel"),
    "conformer_search": ("git", "snakemake", "xtb", "crest"),
    "classical_md": ("git", "snakemake", "packmol", "ambertools", "gromacs", "mpi"),
    "cdft": ("git", "snakemake", "cp2k", "mpi"),
    "embedded_vde": ("git", "snakemake", "orca", "mpi"),
    "plane_wave_gate": ("git", "snakemake", "quantum_espresso", "mpi"),
    "hpc": ("git", "snakemake", "slurm", "mpi"),
    "container": ("git", "apptainer"),
    "production": (
        "git",
        "snakemake",
        "packmol",
        "openbabel",
        "ambertools",
        "gromacs",
        "cp2k",
        "orca",
        "quantum_espresso",
        "xtb",
        "crest",
    ),
}
DOCTOR_REQUIREMENTS: tuple[str, ...] = tuple(REQUIREMENT_ENGINES)


def _first_line(text: str) -> str | None:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped[:300]
    return None


def _find_executable(spec: EngineSpec) -> str | None:
    """Find an engine on PATH or beside the active virtual-environment Python."""
    for name in spec.executables:
        executable = shutil.which(name)
        if executable:
            return executable

    # run.sh deliberately does not mutate PATH when it selects .venv/bin/python.
    # Searching that Python's directory keeps doctor consistent with run.sh.
    python_bin = str(Path(sys.executable).parent)
    for name in spec.executables:
        executable = shutil.which(name, path=python_bin)
        if executable:
            return executable
    return None


def _output_validation_error(spec: EngineSpec, output: str) -> str | None:
    normalized = output.casefold()
    rejected = next(
        (marker for marker in spec.rejected_output_markers if marker.casefold() in normalized),
        None,
    )
    if rejected:
        return f"executable output matched rejected marker {rejected!r}"
    if spec.accepted_output_markers and not any(
        marker.casefold() in normalized for marker in spec.accepted_output_markers
    ):
        return "executable output did not match the expected program signature"
    return None


def required_engine_names(requirements: Sequence[str]) -> frozenset[str]:
    unknown = sorted(set(requirements).difference(REQUIREMENT_ENGINES))
    if unknown:
        available = ", ".join(DOCTOR_REQUIREMENTS)
        raise ValueError(
            f"unknown doctor requirement(s): {', '.join(unknown)}; choose from {available}"
        )
    required: set[str] = set()
    for requirement in requirements:
        required.update(REQUIREMENT_ENGINES[requirement])
    return frozenset(required)


def detect_engine(spec: EngineSpec, timeout_seconds: float | None = None) -> EngineStatus:
    executable = _find_executable(spec)
    if not executable:
        return EngineStatus(
            spec.name,
            False,
            None,
            None,
            spec.category,
            spec.required_for,
            "not found on PATH or beside the active Python interpreter",
        )
    version = None
    error = None
    found = True
    try:
        # Some chemistry executables create CRASH, input_tmp.in, or scratch
        # files even for a version probe. Keep those artifacts out of the Git
        # worktree so a diagnostic cannot block the next safe update.
        with tempfile.TemporaryDirectory(
            prefix=f"solvelec-doctor-{spec.name}-", ignore_cleanup_errors=True
        ) as probe_directory:
            completed = subprocess.run(
                [executable, *spec.version_args],
                cwd=probe_directory,
                capture_output=True,
                text=True,
                timeout=spec.timeout_seconds if timeout_seconds is None else timeout_seconds,
                check=False,
            )
        stdout = completed.stdout or ""
        stderr = completed.stderr or ""
        combined_output = f"{stdout}\n{stderr}"
        version = _first_line(stdout) or _first_line(stderr)
        signature_error = _output_validation_error(spec, combined_output)
        if signature_error:
            found = False
            error = signature_error
        if completed.returncode not in (0, 1) and version is None:
            error = f"version command exited {completed.returncode}"
    except subprocess.TimeoutExpired as exc:
        # A congested compute node or shared filesystem can make an otherwise
        # valid program slow to print its banner.  The executable was still
        # resolved successfully, so preserve that inventory fact and report
        # only that its version/signature remains unconfirmed.
        error = str(exc)
    except (OSError, subprocess.SubprocessError) as exc:
        error = str(exc)
        if spec.accepted_output_markers:
            found = False
    return EngineStatus(
        spec.name, found, executable, version, spec.category, spec.required_for, error
    )


def _doctor_recommendations(
    statuses: Sequence[EngineStatus], missing_required: Sequence[str]
) -> list[str]:
    if not missing_required:
        return []
    by_name = {status.name: status for status in statuses}
    recommendations: list[str] = []
    if "snakemake" in missing_required:
        recommendations.append(
            "Run ./run.sh bootstrap to install Snakemake in the repository .venv."
        )
    orca = by_name.get("orca")
    if "orca" in missing_required:
        if orca and orca.executable and orca.error:
            recommendations.append(
                "The detected 'orca' is not ORCA quantum chemistry (on Linux, /usr/bin/orca is "
                "usually the desktop screen reader). Load the ORCA module or prepend the real "
                "ORCA installation directory to PATH."
            )
        else:
            recommendations.append(
                "Load the licensed ORCA quantum-chemistry module and ensure its installation "
                "directory is on PATH."
            )
    chemistry_missing = [
        name
        for name in missing_required
        if by_name.get(name) and by_name[name].category == "chemistry" and name != "orca"
    ]
    if chemistry_missing:
        recommendations.append(
            "Load the site modules for the missing chemistry engines, or build the bundled "
            "Apptainer environment for redistributable open-source engines."
        )
    return recommendations


def doctor_report(
    specs: Sequence[EngineSpec] = ENGINE_SPECS,
    requirements: Sequence[str] = (),
) -> dict[str, object]:
    selected_requirements = tuple(dict.fromkeys(requirements))
    required_names = required_engine_names(selected_requirements)
    statuses = [detect_engine(spec) for spec in specs]
    engine_records: list[dict[str, object]] = []
    for status in statuses:
        record: dict[str, object] = dict(status.as_dict())
        record["required"] = status.name in required_names
        engine_records.append(record)
    missing_required = sorted(
        status.name for status in statuses if status.name in required_names and not status.found
    )
    return {
        "python": {"version": sys.version.split()[0], "executable": sys.executable},
        "platform": platform.platform(),
        "requirements": list(selected_requirements),
        "ready": not missing_required,
        "missing_required": missing_required,
        "engines": engine_records,
        "missing": [status.name for status in statuses if not status.found],
        "recommendations": _doctor_recommendations(statuses, missing_required),
    }
