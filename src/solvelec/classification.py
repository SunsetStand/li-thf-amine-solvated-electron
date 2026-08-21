from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass


@dataclass(frozen=True)
class LocalizationMetrics:
    electron_count: float
    li_spin: float
    max_atomic_spin: float
    interstitial_fraction: float


@dataclass(frozen=True)
class LocalizationResult:
    label: str
    confidence: str
    reasons: tuple[str, ...]

    def as_dict(self) -> dict[str, str | list[str]]:
        return {
            "label": self.label,
            "confidence": self.confidence,
            "reasons": list(self.reasons),
        }


def classify_localization(
    metrics: LocalizationMetrics, thresholds: Mapping[str, float]
) -> LocalizationResult:
    lower = float(thresholds["electron_count_min"])
    upper = float(thresholds["electron_count_max"])
    li_limit = float(thresholds["li_spin_collapse"])
    atom_limit = float(thresholds["atomic_spin_localized"])
    cavity_limit = float(thresholds["cavity_interstitial_fraction"])

    reasons: list[str] = []
    if not lower <= metrics.electron_count <= upper:
        reasons.append(
            f"integrated spin {metrics.electron_count:.3f} is outside [{lower:.3f}, {upper:.3f}]"
        )
        return LocalizationResult("invalid_spin_integral", "high", tuple(reasons))
    if abs(metrics.li_spin) >= li_limit:
        reasons.append(f"Li spin {metrics.li_spin:.3f} exceeds collapse threshold {li_limit:.3f}")
        return LocalizationResult("li_atomic_or_contact", "high", tuple(reasons))
    if metrics.max_atomic_spin >= atom_limit:
        reasons.append(
            f"maximum atomic spin {metrics.max_atomic_spin:.3f} exceeds {atom_limit:.3f}"
        )
        return LocalizationResult("molecular_anion", "high", tuple(reasons))
    if metrics.interstitial_fraction >= cavity_limit:
        reasons.append(
            "interstitial spin fraction "
            f"{metrics.interstitial_fraction:.3f} exceeds {cavity_limit:.3f}"
        )
        return LocalizationResult("cavity_electron", "medium", tuple(reasons))
    reasons.append("no calibrated localization class passes all thresholds")
    return LocalizationResult("ambiguous", "low", tuple(reasons))
