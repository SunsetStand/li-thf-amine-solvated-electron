# Classical solvent pilot

`classical_pilot` generates the first statistically checked liquid ensembles.
It is deliberately restricted to the `pilot` campaign: pure THF and THF/EDA
at the 1.5 M seed composition, each with three independent replicas. The
trajectories contain neutral solvents only. They do not contain Li or an
explicit excess electron and therefore cannot establish an electron cavity.

## Protocol

Each replica uses GAFF2/AM1-BCC parameters and follows:

1. energy minimization;
2. 0.5 ns NVT at 298.15 K;
3. 5 ns NPT at 298.15 K and 1 bar;
4. 20 ns NPT production.

The time step is 2 fs with hydrogen-bond constraints. PME treats electrostatics,
the 0.9 nm Lennard-Jones cutoff includes a long-range dispersion correction,
and stochastic velocity rescaling plus stochastic cell rescaling sample the
temperature and pressure ensembles. GROMACS 2023 documents `C-rescale` as
appropriate for both equilibration and production and documents checkpoint
continuation with `mdrun -cpi`:

- <https://manual.gromacs.org/documentation/2023.1/user-guide/mdp-options.html>
- <https://manual.gromacs.org/documentation/2023.2/user-guide/managing-simulations.html>

Velocity seeds are deterministic functions of system ID and replica. Different
systems and replicas therefore receive distinct seeds, while a repeated run
from the same commit remains reproducible.

## Required order on TMC-AMD

First validate the two-solvent parameter and topology path:

```bash
./run.sh dry-run --campaign mixed_smoke --target classical_smoke
./run.sh submit --campaign mixed_smoke --target classical_smoke
```

After that target succeeds, preview and submit the six-replica pilot:

```bash
./run.sh dry-run --campaign pilot --target classical_pilot
./run.sh submit --campaign pilot --target classical_pilot
./run.sh queue
```

All commands are Slurm submissions on TMC. The six replica chains can run in
parallel, but EM, NVT, NPT, and production remain ordered within each replica.
Every child requests four CPUs, following the site allocation policy.

## Safe continuation

Long MD stages write partial engine files below a stage-local `.resume/`
directory. A new Slurm allocation checks hashes of the MDP, coordinates,
topology, incoming checkpoint, thread count, and checkpoint interval. It uses
`-cpi` and `-append` only when all hashes match. Changed inputs discard the
stale stage-local resume state and begin that stage again.

Do not delete `.resume/` after a timeout or cancellation. Submit:

```bash
./run.sh resume --campaign pilot --target classical_pilot
```

Successful stages atomically promote their files out of `.resume/`; Snakemake
then treats them as normal completed outputs.

## Acceptance report

Each production trajectory records coordinates every 10 ps. Its validation
JSON checks:

- the GROMACS normal-completion marker;
- at least 98% of the requested 20 ns trajectory span;
- no more than 2% relative density difference between trajectory halves;
- achieved amine concentration within 0.05 M of the target.

The campaign summary additionally requires the three replica mean densities
for each system to span no more than 3% of their mean. Reports are retained even
when a gate fails:

```text
pilot/classical/<system>/r<replica>/pilot/validation.json
pilot/classical_pilot.validation.json
```

`pilot/classical_pilot.done` is created only when every check passes. A failed
concentration check reports a bounded suggested amine count for a subsequent
NPT count-refinement iteration; it never silently changes the composition.

## Scientific boundary

Passing this target establishes solvent density and composition sampling for
pure THF and THF/EDA. It does not pass the Li force-field gate or any electronic
structure gate. The accepted frames become candidates for Li insertion and
CP2K cDFT benchmarks in the next workflow phase.

## Amine-series charge-rounding repair

The `amine_series_pilot` trajectories also sample **neutral THF and neutral
amines only**. The classical force field has neither Li nor an excess electron.
For 1,2-PDA, 1,3-PDA, and TMEDA, the original AM1-BCC MOL2 atom charges can sum
to about 0.001 elementary charge per molecule. Replication across a periodic
box then triggers GROMACS's net-charge/PME warning before minimization starts.
The topology builder keeps the original MOL2 files untouched, writes a local
charge-balanced copy for TLeap only when the discrepancy is at most 0.002 e,
and rejects larger discrepancies. It records each correction in the build
manifest and checks the TLeap residue charges before conversion to GROMACS.

To advance only the six affected systems, first dry-run and submit the EM gate.
Inspect the resulting topology charges and GROMACS logs before submitting their
pilot chains. Both targets use the normal Slurm controller and profile:

```bash
./run.sh dry-run --campaign amine_series_pilot --target amine_series_charge_repair_em
./run.sh submit --campaign amine_series_pilot --target amine_series_charge_repair_em
./run.sh dry-run --campaign amine_series_pilot --target amine_series_charge_repair
./run.sh submit --campaign amine_series_pilot --target amine_series_charge_repair
```

Previous `*.failed` diagnostics are retained even after a successful retry;
later failures are stored under `failure_history/` without overwriting them.

## Two-system 3 M concentration refinement

The first 1,2-PDA and TMEDA 3 M pilots used 21 and 28 amines with 64 THF.
Their measured production concentrations were 3.1018 and 3.1147 M, exceeding
the unchanged 3.00 ± 0.05 M acceptance window. Their validation records remain
failed and are not promoted to analysis or electronic snapshots.

The separate `amine_series_count_refinement` campaign tries 20 1,2-PDA and
27 TMEDA molecules, respectively. It retains 64 THF, the same target and
validation tolerance, and stores new runs under its own campaign directory.
Both the Packmol input and TLeap topology use the overridden counts. The
original `amine_series_pilot` artifacts are preserved for audit.

Preview the two-system DAG and submit through the normal Slurm controller:

```bash
./run.sh dry-run --campaign amine_series_count_refinement --target classical_pilot
./run.sh submit --campaign amine_series_count_refinement --target classical_pilot
```

Inspect the completed job and both replica validation records before using
either trajectory downstream. In particular, NPT may shrink the smaller boxes,
so the new achieved concentrations cannot be guaranteed from the old volume.
If one still misses ±0.05 M, retain its failed result and reconsider box size
before further refinement. Stage B3 source selection remains unchanged until
the new campaign passes validation.
