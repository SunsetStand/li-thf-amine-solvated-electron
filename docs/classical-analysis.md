# Stage A: trajectory analysis and solvent snapshot bank

Stage A turns the completed solvent-only `classical_pilot` trajectories into a
small, auditable descriptor set and a deterministic snapshot bank. It does not
insert lithium, create an excess electron, or run an electronic-structure
engine. Those are later, independently gated stages.

## Scope and execution order

The target is intentionally limited to the `pilot` campaign: pure THF and
THF/EDA 1.5 M, three replicas each. On TMC-AMD, both commands below submit a
controller to Slurm; all analysis and export work is performed by Slurm child
jobs.

```bash
./run.sh dry-run --campaign pilot --target classical_analysis
./run.sh submit --campaign pilot --target classical_analysis
```

Wait for the first controller and all child jobs to finish. Only after
`classical_analysis.summary.json` reports `"ready": true` should the snapshot
bank be submitted:

```bash
./run.sh dry-run --campaign pilot --target snapshot_bank
./run.sh submit --campaign pilot --target snapshot_bank
```

The dependency graph also enforces this order. A failed analysis gate prevents
all snapshot-export rules from starting.

The completed pilot is treated as an immutable upstream data product. Its
TPR/XTC/validation paths are checked and checksummed at analysis runtime but are
not producer dependencies of the Stage-A DAG. This boundary prevents a newer
analysis setting or source file from rescheduling the already validated 20 ns
MD. If a pilot trajectory is intentionally replaced, start a new campaign or
explicitly invalidate its Stage-A outputs before analysis.

## Analysis performed

Every replica is sampled at 100 ps spacing over the full 20 ns production
trajectory. The analysis records:

- instantaneous density and periodic cell volume;
- THF O--O RDFs;
- for the mixture, EDA N--THF O and intermolecular EDA N--N RDFs;
- EDA N--THF O contact counts at 3.5 angstrom;
- EDA N--H...O hydrogen bonds using a 3.5 angstrom donor--acceptor cutoff and
  150 degree angular cutoff;
- the largest periodic heavy-atom van-der-Waals surface-clearance proxy on a
  deterministic grid with local refinement;
- integrated autocorrelation times and effective sample sizes for the scalar
  descriptors used by the readiness gate.

All distances use the current periodic cell. The free-volume value is only a
geometric seed-ranking proxy. It is not electron density, spin density, a
binding energy, or evidence that a solvated electron exists.

Each replica must retain at least 150 analysis frames, have at least five
effective samples for every required scalar descriptor, preserve a positive
void proxy, and depend on a passing classical-pilot validation record. Failure
is retained as a machine-readable JSON record instead of being silently
discarded.

## Deterministic snapshot selection

Stage A exports exactly one representative snapshot from each replica. The
candidate pool is restricted to the second half of production (elapsed time at
least 10 ns). A robust median/MAD-normalized medoid is selected from density
and void radius; the mixed solvent additionally uses contact and hydrogen-bond
counts. Selection is deterministic for unchanged inputs and configuration.

The exported coordinates are wrapped by whole residues where topology data
allow it. Every metadata record states `li_atom_present: false` and preserves
SHA-256 checksums for the source trajectory, topology, analysis, XYZ, and cell
include. No inference about Li coordination or electron localization is valid
from these structures alone.

## Output layout and checks

Large outputs remain below the configured storage run root, normally:

```text
/data/home/storage/Backup_Data/$USER/li-thf-amine-solvated-electron/runs/pilot/
```

The small Stage-A products are:

```text
analysis/<system>/r<replica>/analysis.json
analysis/<system>/r<replica>/timeseries.csv
analysis/<system>/r<replica>/rdf.csv
analysis/<system>/r<replica>/snapshot/representative.xyz
analysis/<system>/r<replica>/snapshot/representative.cell.inc
analysis/<system>/r<replica>/snapshot/metadata.json
classical_analysis.summary.json
classical_analysis.done
snapshot_bank.summary.json
snapshot_bank.done
```

A successful first target has six ready analysis records, two systems with
replicas `[1, 2, 3]`, and `ready: true` in
`classical_analysis.summary.json`. A successful second target has six ready,
unique-replica metadata records and `ready: true` in
`snapshot_bank.summary.json`. The `.done` files contain the SHA-256 hash of the
summary that passed the gate.

Use the normal scheduler and log helpers while jobs run:

```bash
./run.sh queue
./run.sh logs submit
./run.sh status --campaign pilot --target classical_analysis
./run.sh status --campaign pilot --target snapshot_bank
```

Do not invoke the Python analysis script or Snakemake directly on the TMC login
node.
