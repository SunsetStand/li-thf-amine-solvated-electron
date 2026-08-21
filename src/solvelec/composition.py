from __future__ import annotations

import math

AVOGADRO_MOL_INV = 6.022_140_76e23
NM3_TO_L = 1.0e-24


def _round_half_up(value: float) -> int:
    return int(math.floor(value + 0.5))


def concentration_molar(molecule_count: int, volume_nm3: float) -> float:
    if molecule_count < 0:
        raise ValueError("molecule_count must be non-negative")
    if volume_nm3 <= 0:
        raise ValueError("volume_nm3 must be positive")
    return molecule_count / (AVOGADRO_MOL_INV * volume_nm3 * NM3_TO_L)


def count_for_volume(target_concentration_m: float, volume_nm3: float) -> int:
    if target_concentration_m < 0:
        raise ValueError("target_concentration_m must be non-negative")
    if volume_nm3 <= 0:
        raise ValueError("volume_nm3 must be positive")
    exact = target_concentration_m * AVOGADRO_MOL_INV * volume_nm3 * NM3_TO_L
    return _round_half_up(exact)


def initial_amine_count(
    target_concentration_m: float,
    thf_count: int,
    thf_molar_volume_l_mol: float,
    amine_molar_volume_l_mol: float,
) -> int:
    """Estimate count while accounting for the additive's own molar volume.

    The estimate assumes ideal additive volumes and is only a seed for the
    measured-volume NPT iteration.
    """

    if target_concentration_m == 0:
        return 0
    if target_concentration_m < 0 or thf_count <= 0:
        raise ValueError("target concentration and THF count must be non-negative/positive")
    if thf_molar_volume_l_mol <= 0 or amine_molar_volume_l_mol <= 0:
        raise ValueError("molar volumes must be positive")
    denominator = 1.0 - target_concentration_m * amine_molar_volume_l_mol
    if denominator <= 0:
        raise ValueError("target concentration is incompatible with the seed molar volume")
    exact = target_concentration_m * thf_count * thf_molar_volume_l_mol / denominator
    return max(1, _round_half_up(exact))


def suggest_count_after_npt(
    target_concentration_m: float,
    measured_volume_nm3: float,
    current_count: int,
    max_change: int = 3,
) -> int:
    """Suggest a bounded discrete update after measuring an NPT volume."""

    if current_count < 0 or max_change < 1:
        raise ValueError("current_count must be non-negative and max_change positive")
    desired = count_for_volume(target_concentration_m, measured_volume_nm3)
    delta = max(-max_change, min(max_change, desired - current_count))
    return current_count + delta


def mole_fraction(amine_count: int, thf_count: int) -> float:
    if amine_count < 0 or thf_count < 0 or amine_count + thf_count == 0:
        raise ValueError("counts must be non-negative and not both zero")
    return amine_count / (amine_count + thf_count)


def enrichment_factor(local_amine_fraction: float, bulk_amine_fraction: float) -> float:
    if not 0 <= local_amine_fraction <= 1:
        raise ValueError("local fraction must be in [0, 1]")
    if not 0 < bulk_amine_fraction <= 1:
        raise ValueError("bulk fraction must be in (0, 1]")
    return local_amine_fraction / bulk_amine_fraction


def log_odds_enrichment(local_amine_fraction: float, bulk_amine_fraction: float) -> float:
    if not 0 < local_amine_fraction < 1 or not 0 < bulk_amine_fraction < 1:
        raise ValueError("both fractions must be strictly between zero and one")
    local_odds = local_amine_fraction / (1 - local_amine_fraction)
    bulk_odds = bulk_amine_fraction / (1 - bulk_amine_fraction)
    return math.log(local_odds / bulk_odds)
