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

## Stage B2C: pre-existing preferential solvation

Stage B2C tests a narrower, physically cleaner follow-up than fixing an
electron during geometry optimization: do ordinary solvent fluctuations that
are already locally EDA-rich accept an unconstrained excess electron more
favorably than EDA-poor fluctuations? It compares three solvent ensembles:

- pure THF, as a composition-degenerate reference;
- THF/EDA 1.5 M;
- THF/EDA 3 M.

The accepted `pilot` bank remains unchanged. A separate `eda3m_pilot`
campaign generates the 3 M classical ensemble, Stage-A snapshot, and Stage-B
void bank. Stage B2C then reads the three accepted banks as immutable inputs
and publishes its combined smoke result below `runs/pilot/`.

For each system and replica-1 snapshot, all ranked voids are characterized by
the mole fraction of EDA molecular heavy-atom centers within 4, 6, and 8
angstrom. The 6-angstrom shell deterministically selects one `eda_rich` and
one `eda_poor` basis seed. Pure THF uses two distinct voids whose EDA fraction
is identically zero; it cannot express an EDA-rich preference.

Each seed receives two CP2K `RUN_TYPE ENERGY` calculations on identical fixed
nuclei and the same ghost basis center:

- a neutral closed-shell singlet;
- an unconstrained charge -1 doublet, with electron and spin-density cubes.

The within-seed proxy `E(neutral) - E(anion)` is reported in eV. A larger value
means more favorable vertical electron attachment within that matched pair.
It is not valid to compare the raw total energies of different compositions.
The analysis also recomputes local EDA content at the final positive-spin
centroid and sums positive spin fractions over EDA and THF molecules.

### Slurm sequence

First create and accept the missing 3 M EDA source bank:

```bash
./run.sh dry-run --campaign eda3m_pilot --target classical_pilot
./run.sh submit --campaign eda3m_pilot --target classical_pilot
./run.sh submit --campaign eda3m_pilot --target classical_analysis
./run.sh submit --campaign eda3m_pilot --target snapshot_bank
./run.sh submit --campaign eda3m_pilot --target stage_b_candidates
```

After each controller completes and its summary is `ready: true`, preview and
submit the combined B2C smoke:

```bash
./run.sh dry-run --campaign pilot --target stage_b2c_preferential_smoke
./run.sh inspect
./run.sh submit --campaign pilot --target stage_b2c_preferential_smoke
```

With all three candidate banks accepted, the incremental B2C dry-run contains
29 jobs: 3 seed selectors, 6 paired renderers, 12 CP2K single points, 6 pair
analyses, 1 summary, and 1 checksum gate. It must contain exactly six neutral
and six anion CP2K jobs and no GROMACS job.

Success requires `stage_b2c_preferential_smoke.summary.json` with `ready: true`
and a matching `stage_b2c_preferential_smoke.done`. This remains a numerical
smoke screen. PBE delocalization error, the uncorrected charged periodic cell,
one snapshot per composition, and the absence of Li/PFAS prohibit a production
stability, ionization, or kinetic conclusion. It tests pre-existing solvent
selection only; electron-induced nuclear reorganization is deliberately not
claimed.
