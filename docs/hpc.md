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
./run.sh doctor --strict-engines
./run.sh test
./run.sh dry-run --campaign pilot
./run.sh submit --campaign pilot --profile private-mygroup
```

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
