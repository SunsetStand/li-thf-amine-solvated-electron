# Adding systems, concentrations, and methods

## New amine

Add one entry to `configs/systems.yaml` using JSON syntax (JSON is valid YAML):

```json
"new": {
  "display_name": "full chemical name",
  "smiles": "canonical SMILES",
  "residue": "NEW",
  "molar_mass_g_mol": 100.0,
  "molar_volume_l_mol": 0.100,
  "charge": 0
}
```

Then add system IDs such as `new_1p5m` to a campaign. The molar volume is only
an initial packing estimate. The achieved concentration must be recomputed from
the NPT volume.

No Stage-B code contains an amine name. Once the new campaign has passed its
classical and snapshot-bank gates, the same configured void search, Li/cavity
candidate classes, CP2K renderer, and summary logic are reused. The molecule's
SMILES, residue, charge, mass, and volume remain the only chemistry-specific
catalog entry.

## New concentration

System IDs use `p` as the decimal separator: 0.75 M becomes `new_0p75m`. The
parser accepts arbitrary non-negative decimal molarities.

## One-step execution boundary

For an accepted snapshot bank, the portable Stage-B entry point is always:

```bash
./run.sh submit --campaign <campaign> --target stage_b
```

`configs/methods.yaml` controls the candidate distances, number of searched
void sites, selected smoke candidate, and smoke replicas per system. Changing
amine identity or concentration does not require editing a rule or Python
module. Earlier liquid preparation is intentionally kept as a separate gated
handoff until the electronic method passes the pure-THF/EDA validation ladder;
only then should an umbrella target join the complete liquid-to-VDE workflow.

## New method

Add a named block in `configs/methods.yaml`, a renderer/template, and tests that
assert charge, multiplicity, basis, and convergence behavior. A new engine must
never be a silent fallback for a missing requested engine.

## New regression fixture

Small text outputs and tiny grids belong in `tests/fixtures`. Do not commit
licensed program binaries or production data. Golden values must state whether
they test parsing/numerics or a calibrated scientific result.
