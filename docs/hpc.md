# HPC setup

## Supported execution modes

1. Host modules plus a local Python environment.
2. Pixi environment containing open-source engines.
3. Apptainer image built from `containers/Apptainer.def`.

ORCA and VASP are never redistributed. Point the workflow at site-provided
executables/modules.

On an HPC login node with `sbatch` available, `./run.sh bootstrap` installs the
`slurm` optional extra automatically. For an already-created environment use
`python -m pip install -e '.[workflow,slurm]'`.

## SLURM

Copy the full profile without tracking site details:

```bash
cp -r configs/profiles/slurm configs/profiles/private-mygroup
```

Edit `configs/profiles/private-mygroup/config.v9+.yaml` to add the executor
plugin's `slurm_account` and `slurm_partition` keys. Directories matching
`private-*` are gitignored, so credentials and site details do not drift into
GitHub. Load site modules for CP2K, GROMACS, ORCA, and MPI before submitting.
Snakemake's executor plugin maps rule resources to SLURM options; see the
[plugin documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).

Run:

```bash
./run.sh doctor --require input_bundle --require hpc
./run.sh test
./run.sh dry-run --campaign pilot
./run.sh submit --campaign pilot --profile private-mygroup
```

Then validate each expensive stage after loading its site modules:

```bash
./run.sh doctor --require classical_md
./run.sh doctor --require cdft
./run.sh doctor --require embedded_vde
./run.sh doctor --require plane_wave_gate
```

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
