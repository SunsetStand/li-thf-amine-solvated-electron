# Stage B3: amine-series and N-H association screen

Stage B3 asks a narrower question than lithium ionization thermodynamics:

> In fixed THF-rich solvent fluctuations, do N-H-facing local environments and
> amine identity change the vertical stabilization and spatial distribution of
> an unconstrained excess electron?

The workflow remains Li-free and PFAS-free. It must not be interpreted as a
free energy, a rate constant, or the thermodynamics of
`Li(s) -> Li+(solv) + e-(solv)`.

## Panel and controls

The bulk solvent is always THF. The environment panel contains pure THF and
1.5/3 M EDA, 1,2-PDA, 1,3-PDA, DETA, and TMEDA. TMEDA has no N-H bonds and is
the explicit negative control. Covalent N-H pairs are inferred from each
snapshot; the configured molecular counts (EDA 4, 1,2-PDA 4, 1,3-PDA 4, DETA
5, TMEDA 0) are audit expectations rather than atom-name shortcuts.

An N-H group is called "inward-facing" when the periodic N-H-center angle is at
least 120 degrees. Distances, angles, shell counts, and a smooth ranking score
are reported. These are geometric electron-proton association descriptors,
not proof of a classical hydrogen bond.

## Staged targets

### 1. Reanalyze accepted B2C cubes (no new CP2K)

```bash
./run.sh submit --campaign pilot --target stage_b3_nh_reanalysis
```

This adds N-H descriptors to the existing pure-THF/EDA B2C ensemble. It is the
fastest first result and cannot silently reschedule the accepted CP2K jobs.

### 2. Prepare missing solvent ensembles

```bash
./run.sh submit --campaign amine_series_pilot --target classical_pilot
./run.sh submit --campaign amine_series_pilot --target classical_analysis
./run.sh submit --campaign amine_series_pilot --target snapshot_bank
./run.sh submit --campaign amine_series_pilot --target stage_b_candidates
```

The dedicated source campaign contains one replica for each of the eight
previously missing 1,2-PDA, 1,3-PDA, DETA, and TMEDA concentration systems.

### 3. Build the eleven-system environment panel

```bash
./run.sh submit --campaign pilot --target stage_b3_amine_environment
```

This target is lightweight. It selects separated `nh_facing` and `nh_control`
sites from existing voids and reports their N-H geometry and local composition.
For pure THF and TMEDA the N-H contrast is explicitly marked degenerate.

### 4. Run the minimal electronic pilot

```bash
./run.sh submit --campaign pilot --target stage_b3_amine_pilot
```

The pilot contains DETA 3 M and TMEDA 1.5/3 M, two sites per system, and one
neutral/anion fixed-geometry pair per site: 12 CP2K single points total. CP2K
jobs are serialized and request 12 ranks, keeping the controller plus child
budget at 16 CPUs or less.

The main paired quantity is

`[E(neutral)-E(anion)]_nh_facing - [E(neutral)-E(anion)]_nh_control`.

A larger positive value means that the N-H-facing fluctuation has the larger
uncorrected vertical attachment proxy. The result remains exploratory because
PBE can over-delocalize anions and the charged periodic cells are not
finite-size corrected.

### 5. Combined audit gate

```bash
./run.sh submit --campaign pilot --target stage_b3_amine_series
```

The combined gate requires the old-cube reanalysis, all-system environment
panel, and DETA/TMEDA electronic pilot. Only after this gate should a very small
representative subset be promoted to a hybrid/range-separated and finite-size
sensitivity check.

## Primary outputs

- `stage_b3_nh_reanalysis.summary.json`
- `stage_b3_amine_environment.summary.json`
- `stage_b3_amine_pilot.summary.json`
- `stage_b3_amine_series.summary.json`
- matching checksum `.done` gates

If a GROMACS stage fails, its engine logs are copied to undeclared
`*.failed` files beside that stage, and `failure.json` records the Slurm job,
commands, retained paths, and error. These diagnostics survive Snakemake's
automatic removal of incomplete declared outputs.

All large numerical data remain under the configured storage-backed run root.
