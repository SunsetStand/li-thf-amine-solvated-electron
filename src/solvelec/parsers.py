from __future__ import annotations

import re
from dataclasses import dataclass
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


_CP2K_ENERGY = re.compile(
    r"ENERGY\|\s+Total FORCE_EVAL.*?energy.*?:\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)"
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


def parse_output(engine: str, path: str | Path) -> EngineResult:
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    if engine == "cp2k":
        return parse_cp2k_text(text)
    if engine == "orca":
        return parse_orca_text(text)
    raise ValueError(f"Unsupported engine parser {engine!r}")


def vertical_detachment_energy_ev(neutral_doublet_hartree: float, detached_hartree: float) -> float:
    return (detached_hartree - neutral_doublet_hartree) * HARTREE_TO_EV
