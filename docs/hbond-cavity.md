# Neat en and 1:1 THF:en cavity hydrogen-bond workflow

This Stage-A extension compares two solvent-only liquids at 298.15 K and 1 bar:

- `pure_eda`: 64 ethylenediamine molecules and no THF;
- `eda_1to1`: 32 THF plus 32 ethylenediamine molecules, exactly 1:1 by
  molecule count and therefore by molar ratio.

Keeping 64 total solvent molecules makes the structural comparison more direct
than adding 64 EDA molecules to the existing 64-THF box. The initial ideal-volume
box lengths are approximately 19.23 and 19.89 angstrom, respectively; NPT
determines the sampled density and cell volume.

## One-step execution

On TMC-AMD these commands only submit Slurm controllers. Every chemistry and
trajectory-analysis rule runs as a Slurm child job.

```bash
./run.sh dry-run --campaign hbond_smoke --target classical_smoke
./run.sh submit --campaign hbond_smoke --target classical_smoke

./run.sh dry-run --campaign hbond_pilot --target hbond_stage_a
./run.sh submit --campaign hbond_pilot --target hbond_stage_a
```

The pilot target contains three independent 20 ns replicas per system. It waits
for the complete classical-pilot gate before starting analysis, so the second
submit command is a true end-to-end Stage-A workflow rather than two racing
controllers.

## Hydrogen-bond and cavity definitions

An EDA nitrogen is a donor through its bonded hydrogens. Acceptors are THF
oxygen and nitrogen on a different EDA molecule. A hydrogen bond requires:

- donor--acceptor distance no greater than 3.5 angstrom;
- N--H...acceptor angle at least 150 degrees;
- intramolecular EDA N--H...N pairs are excluded.

The cavity remains the largest heavy-atom van-der-Waals surface-clearance proxy
in each sampled periodic frame. For each accepted hydrogen bond, the analysis
measures the shortest distance from the cavity center to the finite
donor--acceptor segment:

- `cavity_associated`: the segment is within 2.0 angstrom of the cavity surface;
- `cavity_bridging`: the segment intersects the cavity sphere with a 0.25
  angstrom numerical margin.

These are geometric labels for the hydrogen-bond network surrounding an empty
region. They do not assert a chemical bond to the cavity center, an electron
density maximum, or a solvated electron.

## Outputs

Each replica produces:

```text
analysis/<system>/r<replica>/analysis.json
analysis/<system>/r<replica>/timeseries.csv
analysis/<system>/r<replica>/rdf.csv
analysis/<system>/r<replica>/hydrogen_bonds.csv
analysis/<system>/r<replica>/snapshot/representative.xyz
analysis/<system>/r<replica>/snapshot/representative.cell.inc
analysis/<system>/r<replica>/snapshot/cavity_hbonds.json
analysis/<system>/r<replica>/snapshot/cavity_local.pdb
analysis/<system>/r<replica>/snapshot/cavity_hbonds.pml
analysis/<system>/r<replica>/snapshot/cavity_hbonds.svg
analysis/<system>/r<replica>/snapshot/metadata.json
```

The SVG contains XY, XZ, and YZ projections. The transparent cyan circle is the
geometric cavity, orange dashed lines are EDA N--H...N bonds, and teal dashed
lines are EDA N--H...O(THF) bonds. Thick dashed lines satisfy the stricter
bridging definition. The local PDB and PyMOL script provide an interactive 3-D
view using the same audited bond list.

`classical_analysis.summary.json` reports each descriptor's mean, sample
standard deviation across the three replicas, and the three individual replica
means. This includes total EDA--EDA/EDA--THF hydrogen bonds and the
cavity-associated and cavity-bridging subsets.

Campaign-level completion requires:

```text
hbond_pilot/classical_pilot.validation.json  ready=true
hbond_pilot/classical_analysis.summary.json  ready=true
hbond_pilot/snapshot_bank.summary.json       ready=true
hbond_pilot/hbond_stage_a.done               exists
```
