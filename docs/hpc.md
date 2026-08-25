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
Completed jobs disappear from `squeue`; `./run.sh queue` itself is a lightweight
query and creates no output file. Locate persistent logs with:

```bash
./run.sh logs
./run.sh logs bootstrap
./run.sh logs probe
```

Inspect the actual compute-node environment:

```bash
./run.sh probe
./run.sh queue
```

`probe` records the allocation, loaded and available relevant modules,
executable paths, and shared filesystems. It also loads the known CP2K, GROMACS,
and ORCA candidates one at a time in isolated subshells, reporting the resulting
executable paths and dependency modules without mixing their environments.

The validated mappings are stage-specific because their MPI dependencies are
incompatible: CP2K 2023.2 uses OpenMPI 4.1.5, GROMACS 2023 uses 4.1.4, and ORCA
6.1.1 uses 4.1.8. The Slurm driver loads the correct family automatically for
each stage-specific doctor job. A production-wide doctor is rejected on this
host; run the separate gates shown below.

The TMC builds can start slowly on the shared filesystem. Engine probes allow
15 seconds by default and 30 seconds for CP2K and repository Snakemake. ORCA
6.1.1 is recognized both by its normal banner and by its distinctive
`parameterfile` diagnostic when it is intentionally launched without an input;
the GNOME screen reader remains explicitly rejected. If a resolved executable
does not print its banner before the timeout, `doctor` keeps it as found and
records the timeout as an unconfirmed version instead of falsely calling the
program missing.

An optional untracked module file can contain only site-wide modules compatible
with every engine:

```bash
mkdir -p configs/profiles/private-tmc-amd
cp configs/profiles/tmc-amd/modules.sh.example \
  configs/profiles/private-tmc-amd/modules.sh
```

Do not put CP2K, GROMACS, ORCA, or MPI in that private file. In particular,
`/usr/bin/orca` is a desktop screen reader and must not be used as the
quantum-chemistry ORCA.

After configuring modules, submit capability gates:

```bash
./run.sh doctor --require input_bundle --require hpc
./run.sh doctor --require classical_md
./run.sh doctor --require cdft
./run.sh doctor --require embedded_vde
./run.sh doctor --require plane_wave_gate
```

## TMC storage and open-source tools

The approved default project root is
`/data/home/storage/Backup_Data/$USER/li-thf-amine-solvated-electron`. The
storage helper rejects roots outside the current user's `Backup_Data`
subdirectory and creates only these project-owned paths:

```text
runs/
software/conda/envs/
software/manifests/
cache/conda-pkgs/
```

Initialize it through Slurm:

```bash
./run.sh storage-init
./run.sh queue
./run.sh logs storage-init
```

Then submit the open-source installation as a separate Slurm job:

```bash
./run.sh tools-install all
./run.sh queue
./run.sh logs tools-install
```

The installation command loads the public `miniconda3` module only after its
Slurm-side Bash process has started. This prevents Miniconda's `libtinfo` from
being injected into the Bash startup while keeping all Conda work inside the
allocation. It writes exact Conda package manifests under
`software/manifests/`. Re-running the command updates the existing prefixes;
individual recovery targets are `chem`, `ambertools`, and `qe`.
Non-interactive confirmation uses `CONDA_ALWAYS_YES=true` rather than the newer
`--yes` subcommand option, so the workflow remains compatible with the site's
older public Conda installation.

Three environments are deliberate. `chem-tools` provides Packmol, OpenBabel,
xTB, and CREST. `ambertools` selects AmberTools 26.0's Python 3.12 no-MPI build
by its full `cuda_None_nompi_py312*` identifier and is
separate because current packaging constrains coexistence with standalone
Packmol. `qe` contains Quantum ESPRESSO and its own OpenMPI. For a
`classical_md` doctor job, the support-tool paths are added first and the site
GROMACS module is loaded afterward, leaving the site's OpenMPI first on `PATH`.
CP2K and ORCA continue to use only their validated site module families.

After installation, submit one capability stage per job:

```bash
./run.sh doctor --require molecule_generation
./run.sh doctor --require conformer_search
./run.sh doctor --require classical_md
./run.sh doctor --require plane_wave_gate
```

Do not activate these Conda environments manually on the login node.

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
2. Stage-separated Conda environments containing open-source engines.
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
On TMC-AMD it is rejected before execution because it would mix incompatible MPI
families; use one engine-stage requirement per command.

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

On the current TMC account, the repository and command working directory stay
at `/data/home/wangcx/li-thf-amine-solvated-electron`. A future production-data
root under `/data/home/storage` is a separate physical storage location, not a
replacement working directory; the workflow will reference it explicitly.

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
