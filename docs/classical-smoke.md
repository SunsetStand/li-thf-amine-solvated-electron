# Classical smoke chain

`classical_smoke` is the first end-to-end executable target. It verifies the
installed chemistry tools, file conversions, stage-specific module loading,
Slurm child submission, and storage routing. It is not an equilibration or a
production trajectory.

## Scientific boundary

The chain contains 64 neutral THF molecules and no Li. The committed safety
setting `allow_unvalidated_li_forcefield: false` is therefore preserved. The
conventional force field also contains no excess electron. Its output can seed
later electronic-structure work only after the relevant method gates are
implemented and passed.

The smoke schedule is:

1. Open Babel 3D construction from the committed THF SMILES.
2. AmberTools GAFF2 atom typing and AM1-BCC charges; `parmchk2` fills missing
   terms.
3. Packmol construction of the initial cubic liquid box.
4. TLeap assembly followed by ParmEd Amber-to-GROMACS conversion.
5. GROMACS energy minimization, 2 ps NVT, and 2 ps NPT at 298.15 K and 1 bar.

The short MD stages check execution only. Density, concentration, and replica
convergence require the longer protocols in `configs/methods.yaml`.

## Submit on TMC-AMD

After all stage-specific doctor jobs report `ready: true`:

```bash
cd /data/home/wangcx/li-thf-amine-solvated-electron
./run.sh update
./run.sh dry-run --campaign smoke --target classical_smoke
./run.sh submit --campaign smoke --target classical_smoke
./run.sh queue
```

Each command returns immediately after submitting a Slurm controller job. The
controller submits child jobs; no chemistry command is launched on the login
node. The wrapper removes only the parent controller's inherited `SLURM_*`
context before Snakemake starts so its status monitor cannot confuse the parent
with a child; all rules are still submitted through Slurm. Re-run with
`./run.sh resume --campaign smoke --target classical_smoke` to continue missing
or incomplete outputs.

## Outputs and checks

Large outputs are under:

```text
/data/home/storage/Backup_Data/$USER/li-thf-amine-solvated-electron/runs/
  _shared/molecules/thf/
  smoke/specs/pure_thf/r1.json
  smoke/inputs/pure_thf/r1/packmol.inp
  smoke/classical/pure_thf/r1/
    packmol.log
    build/{tleap.log,topol.top,conf.gro,manifest.json}
    em/{em.log,validation.json,manifest.json}
    nvt/{nvt.log,nvt.cpt,validation.json,manifest.json}
    npt/{npt.log,npt.cpt,npt.gro,validation.json,manifest.json}
  smoke/classical_smoke.done
```

Packmol must emit its `Success!` marker, and each GROMACS log must contain
`Finished mdrun` without a fatal marker. Those checks are recorded in the JSON
validation files; process exit status alone is not accepted as success.

Controller logs remain below the repository at `runs/slurm/`. Use:

```bash
./run.sh logs dry-run
./run.sh logs submit
./run.sh status --campaign smoke --target classical_smoke
```

If the chain fails, inspect the first failed child rule in the controller log,
then its stage directory. Do not delete successful upstream outputs; after a
code update, submit the same target and Snakemake will continue from the first
missing or stale product.
