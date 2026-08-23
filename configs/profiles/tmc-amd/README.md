# TMC-AMD profile

This profile routes every Snakemake rule to the `amd` partition. Lightweight
rules reserve four CPUs, following the site's four-CPU allocation guidance.
Future MPI CP2K rules reserve 32 single-threaded ranks.

Site modules are loaded by `configs/slurm/tmc-amd-driver.sbatch`. Put licensed
and version-specific module choices in the ignored file:

```bash
mkdir -p configs/profiles/private-tmc-amd
cp configs/profiles/tmc-amd/modules.sh.example \
  configs/profiles/private-tmc-amd/modules.sh
```

Edit that private file after inspecting the output of `./run.sh probe`.
