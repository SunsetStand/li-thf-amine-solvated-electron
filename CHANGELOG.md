# Changelog

## Unreleased

- Made `doctor` capability checks stage-aware with repeatable `--require` gates.
- Detect Snakemake installed beside the repository virtual-environment Python.
- Reject the Linux desktop screen reader when checking for ORCA quantum chemistry.
- Auto-submit task-like `run.sh` commands when Slurm is available, including the
  Snakemake controller itself.
- Add a TMC-AMD four-CPU profile, compute-node probe, private module hook, and
  hard refusal to start chemistry engines outside Slurm.

## 0.1.0 - 2026-08-21

- Initial executable repository scaffold.
- Added composition, PBC spin-density, localization, parser, rendering, and provenance code.
- Added Snakemake/SLURM workflow skeleton and one-command runner.
- Added unit, integration, and scientific-regression tests.
