# HPC setup

## No-login-node execution policy

The TMC-AMD server forbids running tasks directly when Slurm is available.
Repository commands enforce the following boundary:

- `help`, `queue`, lightweight Git synchronization, and `sbatch` stay on the
  login node.
- `bootstrap`, `probe`, `doctor`, `test`, `dry-run`, `status`, `report`, and
  `matrix` are submitted as short Slurm jobs.
- `submit` and `resume` submit a Slurm-hosted Snakemake controller. The
  controller is allowed to submit child jobs, and every Snakemake rule uses the
  Slurm executor.
- Both the Python CLI and the checked chemistry-engine launcher refuse to start
  outside an allocation whenever `sbatch` is visible or
  `SOLVELEC_REQUIRE_SLURM=1` is set.

On a workstation or in CI, where `sbatch` is absent, the same commands run
locally. This preserves dependency-light development without weakening cluster
safety.

## TMC-AMD first use

The bundled driver follows the supplied site template: partition `amd`, module
paths under `/data/modulefiles`, public `miniconda3`, and CPU allocations in
multiples of four. No account or QOS is assumed.

```bash
git pull --ff-only
./run.sh bootstrap
./run.sh queue
```

The submission message contains the job ID and exact `.out`/`.err` paths under
`runs/slurm/`. Wait for bootstrap to complete before submitting the next check.

Inspect the actual compute-node environment:

```bash
./run.sh probe
./run.sh queue
```

`probe` records the allocation, loaded and available relevant modules,
executable paths, and shared filesystems. Use its output to create the untracked
site module file:

```bash
mkdir -p configs/profiles/private-tmc-amd
cp configs/profiles/tmc-amd/modules.sh.example \
  configs/profiles/private-tmc-amd/modules.sh
```

Edit only module names confirmed by the probe. In particular, `/usr/bin/orca`
is a desktop screen reader and must not be used as the quantum-chemistry ORCA.

After configuring modules, submit capability gates:

```bash
./run.sh doctor --require input_bundle --require hpc
./run.sh doctor --require classical_md
./run.sh doctor --require cdft
./run.sh doctor --require embedded_vde
./run.sh doctor --require plane_wave_gate
```

Finally, submit the current input-generation pilot:

```bash
./run.sh dry-run --campaign pilot
./run.sh submit --campaign pilot
```

Both commands are asynchronous Slurm submissions. The current DAG produces
validated input bundles only; expensive chemistry execution rules are not yet
enabled.

## Supported execution modes

All modes below are launched inside Slurm allocations:

1. Host modules plus the repository `.venv`.
2. Pixi environment containing open-source engines.
3. Apptainer image built from `containers/Apptainer.def`.

ORCA and VASP are never redistributed. Point the workflow at site-provided
executables/modules.

Inside its Slurm allocation, `./run.sh bootstrap` installs the `slurm` optional
extra automatically. Do not invoke the underlying Python installation command
manually on the login node.

## SLURM

For a cluster other than TMC-AMD, copy the full profile without tracking site
details:

```bash
cp -r configs/profiles/slurm configs/profiles/private-mygroup
```

Edit `configs/profiles/private-mygroup/config.v9+.yaml` to add the executor
plugin's `slurm_account` and `slurm_partition` keys. Directories matching
`private-*` are gitignored, so credentials and site details do not drift into
GitHub. Load site modules for CP2K, GROMACS, ORCA, and MPI before submitting.
Snakemake's executor plugin maps rule resources to SLURM options; see the
[plugin documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
The automatic wrapper defaults to the TMC driver; another cluster must also set
`SOLVELEC_SLURM_DRIVER`, `SOLVELEC_SLURM_PARTITION`, and
`SOLVELEC_DEFAULT_PROFILE` to its reviewed site-specific values.

Requirements can be combined in one command. `--require production` (and its
legacy alias `--strict-engines`) checks every chemistry engine, so it is expected
to fail on a host prepared only for input generation or one production stage.

On many Linux systems `/usr/bin/orca` is the GNOME desktop screen reader, not
the ORCA quantum-chemistry program. The environment check rejects that false
positive. Load the licensed ORCA module, or prepend the real ORCA installation
directory to `PATH`, before checking `embedded_vde`.

Snakemake installed by `./run.sh bootstrap` lives inside the repository `.venv`.
The environment check searches beside the active virtual-environment Python, so
it reports the same Snakemake executable that `run.sh` actually uses even when
`.venv/bin` is not globally present on `PATH`.

The initial `input_bundle` target only creates validated specifications and
engine inputs. Expensive chemistry rules should be enabled after G0–G2 and
site-specific engine commands have been reviewed.

## Controller resource overrides

The default TMC-AMD controller reserves four CPUs, 8 GB, and 12 hours. Short
diagnostics reserve four CPUs, 4 GB, and 30 minutes. Site-approved overrides can
be supplied to the lightweight submission wrapper, for example:

```bash
SOLVELEC_CONTROLLER_TIME=24:00:00 ./run.sh submit --campaign pilot
SOLVELEC_SLURM_PARTITION=amd ./run.sh probe
```

These variables modify Slurm requests; they do not permit local execution.

## Storage gate before production

The current input-generation pilot writes only small files below `runs/` in the
repository. Before enabling trajectories, wavefunctions, restart files, or cube
outputs, select a user-owned project directory under `/data/home/storage` and
record its quota/backup policy. Do not guess or automatically create a shared
production location. The `probe` output includes filesystem usage needed for
that decision.

## Restarts and provenance

- Keep GROMACS `.cpt` and CP2K restart/WFN files on scratch or project storage.
- Never commit those files to GitHub.
- Do not regenerate velocities on continuation.
- Every run directory contains a manifest with code commit, dirty state,
  software capabilities, inputs, and SHA256 checksums.
- Pin production campaigns to a release tag; do not run them from a moving `main`.

Different hardware is not expected to produce bitwise-identical trajectories.
Reproducibility means pinned inputs, seeds, software, and statistically defined
acceptance tolerances.
