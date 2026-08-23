# TMC-AMD profile

This profile routes every Snakemake rule to the `amd` partition. Lightweight
rules reserve four CPUs, following the site's four-CPU allocation guidance.
Future MPI CP2K rules reserve 32 single-threaded ranks.

Site modules are loaded by `configs/slurm/tmc-amd-driver.sbatch`. Validated
engine-specific choices are pinned in `configs/slurm/tmc-amd-stage-modules.sh`:

- classical MD: `gromacs/2023` with OpenMPI 4.1.4;
- cDFT: `cp2k/2023.2` with OpenMPI 4.1.5;
- embedded VDE: `orca/6.1.1` with OpenMPI 4.1.8.

Put only additional modules that are compatible with all three engines in the
ignored file:

```bash
mkdir -p configs/profiles/private-tmc-amd
cp configs/profiles/tmc-amd/modules.sh.example \
  configs/profiles/private-tmc-amd/modules.sh
```

Never load an MPI or one of the three pinned engines in that private file.
