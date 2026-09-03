# Stage B: Li/cavity candidates and CP2K numerical smoke

Stage B is the first electron-bearing workflow phase. It turns each accepted,
solvent-only Stage-A snapshot into deterministic Li/cavity starting structures,
then tests the periodic CP2K constrained-DFT execution path on a deliberately
small subset. Passing Stage B is an execution and input-consistency result, not
evidence that a solvated electron exists.

## One command and one optional inspection point

On TMC-AMD, every command below submits a Slurm controller and every rule runs
as a child Slurm job:

```bash
./run.sh doctor --require cdft
./run.sh dry-run --campaign pilot --target stage_b
./run.sh submit --campaign pilot --target stage_b
```

The one-step target first builds and gates all candidates, then launches the
CP2K smoke subset. To inspect candidate geometry before spending CP2K time, run
only:

```bash
./run.sh submit --campaign pilot --target stage_b_candidates
```

The current pilot DAG contains 14 jobs: six candidate builders, one candidate
summary and gate, two CP2K input renderers, two CP2K jobs, one smoke summary,
and one final gate. It must contain no classical-MD, trajectory-analysis, or
snapshot-export rule. Stage-A paths are runtime-validated immutable parameters,
so Stage B cannot reschedule or overwrite the accepted trajectories/snapshots.

## Generic candidate construction

For each periodic snapshot, the workflow ranks non-overlapping local maxima of
the nearest-heavy-atom van-der-Waals surface clearance. It then selects distinct
site pairs closest to the configured Li--cavity distances:

- `compact`: target 3.5 angstrom;
- `separated`: target 6.0 angstrom;
- `distant`: target 9.0 angstrom.

The inserted Li center is atom 1 in every CP2K coordinate file. The second site
is written as one `Gh` ghost basis center. A ghost has basis functions but no
nucleus or charge; it gives the initial interstitial region variational support.
Neither a large geometric void nor density on the ghost basis is by itself an
electron-localization observable. Every candidate records the achieved periodic
distance, both clearance radii, source hashes, cell, atom counts, and the exact
Li/ghost roles.

All search settings and candidate definitions live in `configs/methods.yaml`.
The code reads `system_id` and `amine` from the normal system specification; it
contains no EDA-specific branch. Therefore an accepted snapshot bank for any
catalogued amine+THF campaign uses the same command and code path.

## CP2K smoke scope

Only the `separated` candidate from replica 1 of pure THF and THF/EDA 1.5 M is
run in the current pilot. Each job reserves 32 MPI ranks, 128 GB, and 12 hours.
The method is unrestricted periodic PBE-D3(BJ), a modest MOLOPT basis, a ghost
basis center, and a Becke cDFT constraint with `TARGET 2.0` on Li atom 1.

The numerical smoke uses an explicit cDFT population tolerance of 0.05
electrons, separate from the tighter electronic SCF tolerance. This accepts an
execution-path result only when the Li population is within 2.5% of the
two-electron target. The production PBE0 method has a separate configuration and
is not relaxed by this smoke-only setting.

In CP2K 2023.2, the cDFT charge target is the desired valence-electron
population and core charge is not subtracted. The TMC installation's standard
Li potential is explicitly `GTH-PBE-q3`, so Li+ requires two remaining valence
electrons. The renderer couples that potential to `TARGET 2.0` and rejects an
inconsistent pair. See the
[CP2K 2023.2 cDFT input reference](https://manual.cp2k.org/cp2k-2023_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/QS/CDFT.html)
and [KIND/GHOST reference](https://manual.cp2k.org/cp2k-2023_2-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/KIND.html).

The smoke calculation uses a deliberately loose 0.05-electron outer-cDFT
tolerance. This is limited to 2.5% of the 2.0-electron Li target and is justified
only as a numerical pipeline check. In the archived 0.02-electron attempt, the
pure-THF and EDA/THF trajectories reached pre-pathology population errors of
0.03584 e and 0.04502 e, respectively; pursuit of the tighter threshold then
caused long Newton line searches and severe branch excursions. The smoke gate
independently
parses the final CP2K cDFT target and population, recomputes their absolute
deviation, and requires it to be no greater than the configured 0.05 e. It also
requires SCF convergence, a CP2K energy, and normal termination. The summary
records the measured cDFT population, deviation, tolerance, iteration, and
constraint strength for audit.

For unrestricted calculations, CP2K 2023.2's `E_DENSITY_CUBE` print key writes
both electronic and spin-density cubes; the later localization analysis uses
the spin-density member of that output rather than the unsupported
`SPIN_DENSITY_CUBE` subsection.
Its summary is explicitly labelled
`NUMERICAL_SMOKE_ONLY_NOT_A_LOCALIZATION_RESULT`. PBE self-interaction can
qualitatively change an excess-electron state, and the ghost basis is not yet
converged. No VDE, localization class, stability, or amine trend may be reported
from this smoke calculation.

## Outputs and success checks

Large files are written below the configured storage run root, normally:

```text
/data/home/storage/Backup_Data/$USER/li-thf-amine-solvated-electron/runs/pilot/
```

Important products are:

```text
stage_b/<system>/r<replica>/candidates/manifest.json
stage_b/<system>/r<replica>/candidates/<candidate>/coordinates.xyz
stage_b/<system>/r<replica>/candidates/<candidate>/cell.inc
stage_b/<system>/r<replica>/candidates/<candidate>/metadata.json
stage_b/<system>/r1/smoke/separated/cp2k.inp
stage_b/<system>/r1/smoke/separated/cp2k.out
stage_b_candidates.summary.json
stage_b_cp2k_smoke.summary.json
stage_b.done
```

Success requires `ready: true` in both summaries and a `stage_b.done` whose
hash matches the smoke summary. While jobs run:

```bash
./run.sh queue
./run.sh logs submit
./run.sh status --campaign pilot --target stage_b
```

## What follows

After this smoke gate passes, the next electronic ladder is: basis/cutoff and
ghost-basis calibration on pure THF; PBE0/ADMM fixed-geometry calculations;
spin-density/Hirshfeld localization classification; cDFT release; detached-state
energies and VDE; then replica, functional, and finite-size checks. Only a
candidate that passes those independent gates can enter an amine comparison.
