# TMC-AMD profile

This profile routes every Snakemake rule to the `amd` partition. Lightweight
rules reserve four CPUs, following the site's four-CPU allocation guidance.
Future MPI CP2K rules reserve 32 single-threaded ranks.

The TMC policy also places the Snakemake controller in a Slurm allocation.
`run.sh` strips that parent allocation's `SLURM_*` variables from the Snakemake
subprocess before the Slurm executor starts its status thread. The controller
still uses this profile to submit every rule as a child Slurm job; the stage
guard remains enabled and rejects chemistry outside a child allocation.

The profile explicitly uses `squeue` for child-job status checks. TMC-AMD keeps
completed jobs in the controller for 300 seconds, which the executor plugin
reports as sufficient for reliable `squeue` monitoring. This avoids the site's
incompatible `sacct` account/status query behavior, which can otherwise leave a
controller waiting after successful child jobs have finished.

Site modules are loaded by `configs/slurm/tmc-amd-driver.sbatch`. Validated
engine-specific choices are pinned in `configs/slurm/tmc-amd-stage-modules.sh`:

- classical MD: `gromacs/2023` with OpenMPI 4.1.4;
- cDFT: `cp2k/2023.2` with OpenMPI 4.1.5;
- embedded VDE: `orca/6.1.1` with OpenMPI 4.1.8.

The site CP2K module does not export its data directory. The cDFT stage pins
`CP2K_DATA_DIR=/data/softwares/cp2k/2023.2/data` and verifies the required
basis, potential, and D3 files before starting CP2K.
MPI rules set the standard Snakemake `mpi` resource to `srun`; this exposes the
full Slurm allocation to CP2K instead of nesting `mpirun` inside the executor's
single-task wrapper step.

Put only additional modules that are compatible with all three engines in the
ignored file:

```bash
mkdir -p configs/profiles/private-tmc-amd
cp configs/profiles/tmc-amd/modules.sh.example \
  configs/profiles/private-tmc-amd/modules.sh
```

Never load an MPI or one of the three pinned engines in that private file.
