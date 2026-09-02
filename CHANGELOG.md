# Changelog

## Unreleased

- Made `doctor` capability checks stage-aware with repeatable `--require` gates.
- Detect Snakemake installed beside the repository virtual-environment Python.
- Reject the Linux desktop screen reader when checking for ORCA quantum chemistry.
- Auto-submit task-like `run.sh` commands when Slurm is available, including the
  Snakemake controller itself.
- Add a TMC-AMD four-CPU profile, compute-node probe, private module hook, and
  hard refusal to start chemistry engines outside Slurm.
- Pass the repository root explicitly into Slurm jobs instead of deriving it
  from Slurm's copied `/var/spool/.../slurm_script` path.
- Add a login-safe `logs` command and isolated candidate-module activation to
  the compute-node probe.
- Pin TMC-AMD CP2K, GROMACS, and ORCA modules per workflow stage and reject
  doctor requests that would mix their incompatible OpenMPI versions.
- Recognize the ORCA 6.1.1 no-input signature and lengthen probes for CP2K and
  Snakemake on slow shared filesystems.
- Treat a version-probe timeout as an unconfirmed version rather than a missing
  executable, while still rejecting completed probes with invalid signatures.
- Add Slurm-only TMC storage initialization under each user's `Backup_Data`
  directory and stage-separated Conda installation for chemistry tools,
  no-MPI AmberTools, and Quantum ESPRESSO.
- Match AmberTools 26's full `cuda_None_nompi_py312*` build string instead of
  incorrectly assuming that the build identifier begins with `nompi_`.
- Delay loading Miniconda until after the Slurm-side Bash entry point starts so
  its bundled `libtinfo` is not selected while launching `run.sh`.
- Use the backward-compatible `CONDA_ALWAYS_YES` setting because the TMC Conda
  environment subcommands do not recognize the newer `--yes` option.
- Add Slurm-only Stage-A analysis for the six completed classical-pilot
  trajectories, including periodic structural descriptors, autocorrelation
  gates, and deterministic solvent-only snapshot export.

## 0.1.0 - 2026-08-21

- Initial executable repository scaffold.
- Added composition, PBC spin-density, localization, parser, rendering, and provenance code.
- Added Snakemake/SLURM workflow skeleton and one-command runner.
- Added unit, integration, and scientific-regression tests.
