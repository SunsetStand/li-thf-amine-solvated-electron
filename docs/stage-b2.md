# Stage B2: Li-free intrinsic excess-electron smoke

Stage B2 implements the deliberately Li-free question requested after the
paired Stage-B/B1 analysis: if one electron is added to a frozen solvent
snapshot, without PFAS and without a localization constraint, where does its
unpaired spin settle numerically?

The experimental context recorded with every result is:

- electron generation is reported to be rate-limiting;
- excess electrons accumulate when PFAS is absent;
- PFAS consumes the electrons on a seconds-or-faster timescale;
- PFAS is therefore excluded from this intrinsic reservoir-state smoke test.

These are user-reported experimental observations, not outputs of this
workflow. PFAS capture and reaction kinetics require a later, separate branch.

## Pilot design

The current pilot uses the already accepted Stage-B solvent environments for
`pure_thf` and `eda_1p5m`, replica 1 only. It removes Li and the old ghost
center, preserves the solvent nuclei and periodic cell, and creates two inputs
from the two highest-ranked instantaneous voids:

- `void_01`
- `void_02`

Each input is a charge -1 doublet with one ghost hydrogen basis center at its
seed void. The ghost supplies basis flexibility but has no nucleus or charge.
There is no cDFT block, spatial constraint, PFAS, Li, or nuclear relaxation.
Using two independent seed locations tests whether the result is dominated by
the chosen basis seed.

The four CP2K jobs use periodic unrestricted PBE-D3(BJ), 32 MPI ranks, 128 GB,
and a 12-hour limit. This is intentionally a numerical vertical-injection
smoke test. Semilocal-functional delocalization error, charged-cell finite-size
effects, and missing structural relaxation prevent a stability conclusion.

## Slurm-only execution

On TMC-AMD:

```bash
./run.sh dry-run --campaign pilot --target stage_b2_intrinsic_smoke
./run.sh inspect
./run.sh submit --campaign pilot --target stage_b2_intrinsic_smoke
```

With the accepted Stage-B candidate bank already present, the incremental
dry-run should contain exactly 16 jobs:

- 2 solvent/seed preparation jobs;
- 4 CP2K input renderers;
- 4 CP2K runs;
- 4 cube analyses;
- 1 summary;
- 1 checksum gate.

It must contain no GROMACS job and exactly four
`run_stage_b2_intrinsic_smoke` jobs.

## Acceptance and interpretation

Success requires both files under the configured campaign run root:

```text
stage_b2_intrinsic_smoke.summary.json
stage_b2_intrinsic_smoke.done
```

The summary must have `ready: true`, four ready state records, normal CP2K
termination, converged SCF, signed spin close to one electron, and valid density
cubes. It reports for each seed:

- spin centroid and radius;
- nearest real atom and van-der-Waals surface clearance;
- solvent-molecule, interstitial, and seed-probe spin fractions;
- whether the two initial seed locations reach the same geometric proxy class;
- their final-energy spread as a reproducibility diagnostic.

A passing gate means only that the four Li-free vertical electronic states were
computed and analyzed consistently. It does not prove an equilibrium solvated
electron, a persistent physical cavity, favorable Li ionization, or PFAS
capture kinetics. A production decision needs at least a better
self-interaction-controlled method, cell-size convergence, neutral reference
and detachment/attachment energetics, nuclear relaxation or dynamics, and
independent wavefunction/topological analysis.
