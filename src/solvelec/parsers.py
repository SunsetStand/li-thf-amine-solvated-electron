from __future__ import annotations

import re
from dataclasses import dataclass
from math import isfinite
from pathlib import Path

HARTREE_TO_EV = 27.211_386_245_988


@dataclass(frozen=True)
class EngineResult:
    engine: str
    converged: bool
    normal_termination: bool
    energy_hartree: float | None
    problems: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "engine": self.engine,
            "converged": self.converged,
            "normal_termination": self.normal_termination,
            "energy_hartree": self.energy_hartree,
            "problems": list(self.problems),
        }


@dataclass(frozen=True)
class CdftIteration:
    iteration: int
    rms_gradient: float
    energy_hartree: float
    target_electrons: float
    current_electrons: float
    reported_deviation_electrons: float
    constraint_strength: float


@dataclass(frozen=True)
class CdftConstraintResult:
    converged: bool
    tolerance_electrons: float
    expected_target_electrons: float
    iteration: CdftIteration | None
    deviation_electrons: float | None
    problems: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        iteration = self.iteration
        return {
            "converged": self.converged,
            "tolerance_electrons": self.tolerance_electrons,
            "expected_target_electrons": self.expected_target_electrons,
            "iteration": iteration.iteration if iteration else None,
            "rms_gradient": iteration.rms_gradient if iteration else None,
            "energy_hartree": iteration.energy_hartree if iteration else None,
            "target_electrons": iteration.target_electrons if iteration else None,
            "current_electrons": iteration.current_electrons if iteration else None,
            "deviation_electrons": self.deviation_electrons,
            "reported_deviation_electrons": (
                iteration.reported_deviation_electrons if iteration else None
            ),
            "constraint_strength": iteration.constraint_strength if iteration else None,
            "problems": list(self.problems),
        }


_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
_CP2K_ENERGY = re.compile(
    r"ENERGY\|\s+Total FORCE_EVAL.*?energy.*?:\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)"
)
_CP2K_CDFT_ITERATION = re.compile(
    rf"CDFT SCF iter\s*=\s*(\d+)\s+RMS gradient\s*=\s*({_FLOAT})"
    rf"\s+energy\s*=\s*({_FLOAT}).*?"
    rf"Target value of constraint\s*:\s*({_FLOAT}).*?"
    rf"Current value of constraint\s*:\s*({_FLOAT}).*?"
    rf"Deviation from target\s*:\s*({_FLOAT}).*?"
    rf"Strength of constraint\s*:\s*({_FLOAT})",
    re.DOTALL,
)
_ORCA_ENERGY = re.compile(r"FINAL SINGLE POINT ENERGY\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)")


def parse_cp2k_text(text: str) -> EngineResult:
    problems: list[str] = []
    normal = "PROGRAM ENDED AT" in text
    scf_failure = any(
        marker in text
        for marker in (
            "SCF run NOT converged",
            "SCF run not converged",
            "ABORT",
            "Cholesky decompose failed",
        )
    )
    matches = _CP2K_ENERGY.findall(text)
    energy = float(matches[-1]) if matches else None
    if not normal:
        problems.append("missing CP2K normal-termination marker")
    if scf_failure:
        problems.append("CP2K reported an SCF or fatal failure")
    if energy is None:
        problems.append("no CP2K FORCE_EVAL energy found")
    return EngineResult(
        "cp2k", normal and not scf_failure and energy is not None, normal, energy, tuple(problems)
    )


def parse_cp2k_cdft_iterations(text: str) -> tuple[CdftIteration, ...]:
    """Parse completed CP2K cDFT outer-SCF iterations and their constraint records."""

    return tuple(
        CdftIteration(
            iteration=int(match.group(1)),
            rms_gradient=float(match.group(2)),
            energy_hartree=float(match.group(3)),
            target_electrons=float(match.group(4)),
            current_electrons=float(match.group(5)),
            reported_deviation_electrons=float(match.group(6)),
            constraint_strength=float(match.group(7)),
        )
        for match in _CP2K_CDFT_ITERATION.finditer(text)
    )


def evaluate_cp2k_cdft_constraint(
    text: str, *, expected_target_electrons: float, tolerance_electrons: float
) -> CdftConstraintResult:
    """Independently gate the final cDFT population against a configured tolerance."""

    if not isfinite(expected_target_electrons):
        raise ValueError("expected cDFT target must be finite")
    if not isfinite(tolerance_electrons) or tolerance_electrons <= 0:
        raise ValueError("cDFT tolerance must be finite and positive")
    iterations = parse_cp2k_cdft_iterations(text)
    if not iterations:
        return CdftConstraintResult(
            False,
            tolerance_electrons,
            expected_target_electrons,
            None,
            None,
            ("no completed CP2K cDFT outer-SCF iteration found",),
        )
    final = iterations[-1]
    deviation = final.current_electrons - expected_target_electrons
    problems: list[str] = []
    target_tolerance = 1.0e-10
    if abs(final.target_electrons - expected_target_electrons) > target_tolerance:
        problems.append(
            "CP2K cDFT target "
            f"{final.target_electrons:.12g} differs from configured target "
            f"{expected_target_electrons:.12g}"
        )
    if abs(deviation) > tolerance_electrons:
        problems.append(
            f"final cDFT population error {abs(deviation):.12g} e exceeds "
            f"{tolerance_electrons:.12g} e"
        )
    return CdftConstraintResult(
        not problems,
        tolerance_electrons,
        expected_target_electrons,
        final,
        deviation,
        tuple(problems),
    )


def parse_orca_text(text: str) -> EngineResult:
    problems: list[str] = []
    normal = "ORCA TERMINATED NORMALLY" in text
    failure = any(
        marker in text for marker in ("SCF NOT CONVERGED", "ORCA finished by error termination")
    )
    matches = _ORCA_ENERGY.findall(text)
    energy = float(matches[-1]) if matches else None
    if not normal:
        problems.append("missing ORCA normal-termination marker")
    if failure:
        problems.append("ORCA reported a failure")
    if energy is None:
        problems.append("no ORCA final energy found")
    return EngineResult(
        "orca", normal and not failure and energy is not None, normal, energy, tuple(problems)
    )


def parse_packmol_text(text: str) -> EngineResult:
    normalized = text.casefold()
    normal = "success!" in normalized
    failure = any(marker in normalized for marker in ("error", "stop 171", "segmentation fault"))
    problems: list[str] = []
    if not normal:
        problems.append("missing Packmol success marker")
    if failure:
        problems.append("Packmol reported a fatal failure")
    return EngineResult("packmol", normal and not failure, normal, None, tuple(problems))


def parse_gromacs_text(text: str) -> EngineResult:
    normalized = text.casefold()
    normal = "finished mdrun" in normalized
    failure = any(
        marker in normalized
        for marker in ("fatal error", "segmentation fault", "nan detected", "core dumped")
    )
    problems: list[str] = []
    if not normal:
        problems.append("missing GROMACS Finished mdrun marker")
    if failure:
        problems.append("GROMACS reported a fatal failure")
    return EngineResult("gromacs", normal and not failure, normal, None, tuple(problems))


def parse_output(engine: str, path: str | Path) -> EngineResult:
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    if engine == "cp2k":
        return parse_cp2k_text(text)
    if engine == "orca":
        return parse_orca_text(text)
    if engine == "packmol":
        return parse_packmol_text(text)
    if engine == "gromacs":
        return parse_gromacs_text(text)
    raise ValueError(f"Unsupported engine parser {engine!r}")


def vertical_detachment_energy_ev(neutral_doublet_hartree: float, detached_hartree: float) -> float:
    return (detached_hartree - neutral_doublet_hartree) * HARTREE_TO_EV
