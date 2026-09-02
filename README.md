# Li/THF/amine solvated-electron workflow

Reproducible, safety-gated computational workflow for studying how EDA,
1,2-PDA, 1,3-PDA, DETA, and TMEDA additives change Li⁺ solvation, excess-
electron cavities, Li–electron pairing, and vertical detachment energies in
THF.

The repository deliberately separates three questions that are easy to mix up:

1. Does a calculation contain a genuine interstitial/cavity electron?
2. Is that state stable after the Li⁺ charge constraint is released?
3. How do finite-temperature ensembles change VDE and Li–electron pairing?

Classical MD is used to generate liquid configurations, not to claim that an
electron cavity has already equilibrated. Electron-bearing states are built and
validated with unrestricted hybrid DFT and constrained DFT before they enter any
scientific comparison.

## Quick start

On a workstation, commands run directly when `sbatch` is absent. On the TMC-AMD
server, the same task-like commands automatically submit themselves to Slurm;
only the lightweight shell wrapper and `sbatch` run on the login node:

```bash
git clone https://github.com/SunsetStand/li-thf-amine-solvated-electron.git
cd li-thf-amine-solvated-electron

./run.sh bootstrap
./run.sh queue
```

Wait for bootstrap to finish, then submit compute-node inspection and input-DAG
checks:

```bash
./run.sh probe
./run.sh doctor --require input_bundle
./run.sh dry-run --campaign pilot
```

Each command prints its job ID and log paths under `runs/slurm/`. Do not run the
underlying Python or Snakemake commands manually on the TMC login node.
Use `./run.sh logs probe` or `./run.sh logs bootstrap` to list matching absolute
log paths. `./run.sh queue` is only a live scheduler query and creates no log.

`doctor` reports missing chemistry engines without silently substituting a
different method. A plain `./run.sh doctor` is an inventory; `--require` turns
one or more workflow stages into enforced gates. For example, use
`./run.sh doctor --require classical_md --require hpc` before submitting the
classical-MD stage. `--strict-engines` remains as a compatibility alias for
`--require production`, which intentionally requires every chemistry engine.
On TMC-AMD, run CP2K, GROMACS, and ORCA gates separately because their site
modules provide different OpenMPI versions.

For the validated TMC account, the repository working directory remains
`/data/home/wangcx/li-thf-amine-solvated-electron`. Initialize the separate
large-data and open-source-tool tree on the storage disk, then install the
staged environments inside Slurm:

```bash
./run.sh storage-init
./run.sh tools-install all
```

Wait for each submitted job to finish. The default storage root is
`/data/home/storage/Backup_Data/$USER/li-thf-amine-solvated-electron`; it is not
the command working directory. AmberTools, general chemistry tools, and Quantum
ESPRESSO are kept in separate environments so their package and MPI constraints
cannot replace the site MPI used by GROMACS, CP2K, or ORCA.
`bootstrap` creates or updates `.venv`, installs the declared dependencies, and
runs the complete dependency-light test suite inside its allocation. Later
checks only need `./run.sh test`; `./run.sh update` performs only the lightweight
Git fast-forward locally, then submits bootstrap through Slurm.

To submit the pilot through the bundled SLURM profile:

```bash
./run.sh submit --campaign pilot
```

On TMC-AMD this first submits the Snakemake controller as a four-CPU `amd` job.
That controller then submits every workflow rule as a child Slurm job. Child
rules use `configs/profiles/tmc-amd`; no workflow rule is marked local. Before
Snakemake starts, the wrapper removes the controller allocation's inherited
`SLURM_*` variables from the Snakemake subprocess. Each child receives a fresh
Slurm environment from `sbatch`, while `SOLVELEC_REQUIRE_SLURM=1` remains in
force. This prevents nested-controller status-thread hangs without permitting
login-node or in-allocation local chemistry execution.

For site-specific account/partition settings, copy the SLURM profile to an
untracked `configs/profiles/private-<site>/` directory and select that profile;
see `docs/hpc.md`. Load site chemistry modules before submitting.

The first executable chemistry target is a deliberately short, neutral-solvent
smoke chain:

```bash
./run.sh dry-run --campaign smoke --target classical_smoke
./run.sh submit --campaign smoke --target classical_smoke
```

It parameterizes THF with GAFF2/AM1-BCC, packs 64 molecules, converts the Amber
topology to GROMACS, minimizes, and runs 2 ps NVT plus 2 ps NPT. It proves
software interoperability only; it is not equilibrated scientific sampling.
No Li is included while `allow_unvalidated_li_forcefield` is false. CP2K,
ORCA, and Quantum ESPRESSO execution remain behind their scientific input and
benchmark gates. See `docs/classical-smoke.md`.

After the pure and mixed-solvent smoke chains pass, the explicitly selected
`classical_pilot` target runs pure THF and THF/EDA 1.5 M with three replicas,
using 0.5 ns NVT, 5 ns NPT equilibration, and 20 ns production:

```bash
./run.sh submit --campaign mixed_smoke --target classical_smoke
./run.sh submit --campaign pilot --target classical_pilot
```

The pilot has hash-guarded checkpoint continuation and density, concentration,
trajectory-length, and replica-consistency gates. It remains a solvent-only
ensemble for later electronic-structure seeds. See `docs/classical-pilot.md`.

Stage A analyzes those six completed trajectories without modifying them. It
computes periodic RDF/contact/hydrogen-bond descriptors, estimates descriptor
autocorrelation and effective sample sizes, and measures a heavy-atom
van-der-Waals clearance proxy for interstitial free volume. Run the analysis
gate first, then export one deterministic representative solvent snapshot per
replica only if every analysis record passed:

```bash
./run.sh submit --campaign pilot --target classical_analysis
./run.sh submit --campaign pilot --target snapshot_bank
```

The snapshots contain solvent only: no Li atom and no explicit electron are
added at this stage. See `docs/classical-analysis.md` for output paths,
acceptance checks, and interpretation limits.

## What is implemented in v0.1

- Config-driven 11-system matrix (pure THF plus five amines at 1.5/3.0 M).
- NPT-volume-aware concentration/count calculations.
- Charge/multiplicity state validation for solvated, detached, and collapsed-Li states.
- Periodic Gaussian-cube spin-density integration, centroid, radius, and IPR.
- Cavity/Li/molecular-anion localization classification with explicit uncertainty.
- CP2K PBE0/ADMM/cDFT and ORCA ΔSCF input rendering.
- Packmol execution and checked AmberTools/ParmEd-to-GROMACS topology conversion.
- Neutral-solvent GROMACS smoke and six-replica classical pilot execution with
  checkpoint-safe continuation and quantitative acceptance gates.
- Periodic pilot-trajectory analysis, autocorrelation-aware readiness checks,
  and deterministic solvent-only electronic-structure seed selection.
- Provenance manifests with Git SHA, software versions, inputs, and checksums.
- Dependency-light unit and integration tests using only core package dependencies.

External chemistry binaries and licensed programs are not redistributed. Large
trajectories, wavefunctions, restart files, and cube files stay on the group
filesystem and are referenced by checksummed manifests.

## Scientific gates

Production conclusions require all gates in `docs/methodology.md` to pass. In
particular, a job is not considered scientifically successful merely because
the SCF converged. Spin-density localization, constraint release, basis/plane-
wave agreement, and finite-size behavior are independent acceptance checks.

## Repository map

- `src/solvelec/`: configuration, analysis, engine, and provenance code.
- `configs/`: campaign, systems, methods, and scheduler settings.
- `workflow/`: Snakemake rules, scripts, and input templates.
- `tests/`: unit/integration/scientific regression fixtures.
- `docs/`: methodology, HPC setup, and extension guides.

## License and citation

Code is MIT licensed. Cite the tagged repository release and the primary
method/software papers appropriate to the calculations actually used.
