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

On a Linux workstation or HPC login node:

```bash
git clone https://github.com/SunsetStand/li-thf-amine-solvated-electron.git
cd li-thf-amine-solvated-electron

./run.sh bootstrap
./run.sh doctor
./run.sh dry-run --campaign pilot
```

`doctor` reports missing chemistry engines without silently substituting a
different method. Use `./run.sh doctor --strict-engines` on the production host.
`bootstrap` creates or updates `.venv`, installs the declared dependencies, and
runs the complete dependency-light test suite. Later checks only need
`./run.sh test`; `./run.sh update` fast-forwards the repository, resynchronizes
the environment, and retests it.

To submit the pilot through the bundled SLURM profile:

```bash
./run.sh submit --campaign pilot --profile slurm
```

The bootstrap command installs the SLURM executor plugin only when `sbatch` is
present. This keeps local/Windows dry-runs independent of cluster commands.

For site-specific account/partition settings, copy the SLURM profile to an
untracked `configs/profiles/private-<site>/` directory and select that profile;
see `docs/hpc.md`. Load site chemistry modules before submitting.

## What is implemented in v0.1

- Config-driven 11-system matrix (pure THF plus five amines at 1.5/3.0 M).
- NPT-volume-aware concentration/count calculations.
- Charge/multiplicity state validation for solvated, detached, and collapsed-Li states.
- Periodic Gaussian-cube spin-density integration, centroid, radius, and IPR.
- Cavity/Li/molecular-anion localization classification with explicit uncertainty.
- CP2K PBE0/ADMM/cDFT and ORCA ΔSCF input rendering.
- Packmol input rendering and chemistry-engine capability checks.
- Snakemake input-generation DAG plus local/SLURM profiles.
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
