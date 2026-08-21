from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ElectronicState:
    name: str
    charge: int
    multiplicity: int
    description: str

    @property
    def unpaired_electrons_minimum(self) -> int:
        return self.multiplicity - 1


SOLVATED_ELECTRON = ElectronicState(
    "solvated_electron", 0, 2, "Li+ plus one solvated electron in a neutral doublet cell"
)
DETACHED = ElectronicState(
    "detached", 1, 1, "Same nuclei after vertical removal of the excess electron"
)
SOLVENT_ONLY = ElectronicState("solvent_only", 0, 1, "Closed-shell neutral solvent reference")
ELECTRON_ONLY_CLUSTER = ElectronicState(
    "electron_only_cluster", -1, 2, "Anionic solvent cluster without Li"
)

STATES = {
    state.name: state
    for state in (SOLVATED_ELECTRON, DETACHED, SOLVENT_ONLY, ELECTRON_ONLY_CLUSTER)
}


def get_state(name: str) -> ElectronicState:
    try:
        return STATES[name]
    except KeyError as exc:
        raise ValueError(
            f"Unknown electronic state {name!r}; choose from {sorted(STATES)}"
        ) from exc


def validate_state(name: str, charge: int, multiplicity: int) -> None:
    state = get_state(name)
    if charge != state.charge or multiplicity != state.multiplicity:
        raise ValueError(
            f"State {name} requires charge={state.charge}, multiplicity={state.multiplicity}; "
            f"received charge={charge}, multiplicity={multiplicity}"
        )
