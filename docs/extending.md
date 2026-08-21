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

## New concentration

System IDs use `p` as the decimal separator: 0.75 M becomes `new_0p75m`. The
parser accepts arbitrary non-negative decimal molarities.

## New method

Add a named block in `configs/methods.yaml`, a renderer/template, and tests that
assert charge, multiplicity, basis, and convergence behavior. A new engine must
never be a silent fallback for a missing requested engine.

## New regression fixture

Small text outputs and tiny grids belong in `tests/fixtures`. Do not commit
licensed program binaries or production data. Golden values must state whether
they test parsing/numerics or a calibrated scientific result.
