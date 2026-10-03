# EDA/TMEDA: Li electron transfer and structural relaxation

This calculation compares EDA and TMEDA in THF at approximately 3 M by the
energy required to transfer the Li valence electron into the solvent, including
the energy gained by relaxing the resulting structure. Both electronic states
contain the same atoms and electrons in a neutral, periodic doublet cell.

## States and calculation sequence

| State | Li Becke valence population | Meaning |
|---|---:|---|
| a | 3 | Li⁰ reference diabatic state |
| b | 2 | Li⁺ and an electron in the solvent diabatic state |

Li uses GTH-PBE-q3. The two populations define consistently partitioned cDFT
states. All solvent atoms and Li can relax; the cell vectors remain fixed.
Only the Li population is constrained. The solvent electron distribution is
determined by the electronic calculation.

For each selected solvent snapshot, the workflow:

1. Adds one Li near amine N donors using the same steric-clearance rule for both
   amines. The initial Li–N distance is 2.2 Å. Bidentate positions are tried
   first, with a recorded monodentate fallback when necessary.
2. Optimizes state a to obtain Rₐ.
3. Starts state b from Rₐ and optimizes it to obtain Rᵦ.
4. Evaluates both electronic states at both optimized geometries.

Pair labels put the electronic state first and the geometry second:

| Output pair | Energy |
|---|---|
| aa | Eₐ(Rₐ) |
| ba | Eᵦ(Rₐ) |
| bb | Eᵦ(Rᵦ) |
| ab | Eₐ(Rᵦ) |

The default method is PBE0-D3(BJ)/ADMM with 25% exact exchange and 75% PBE
exchange, TZV2P-MOLOPT-GTH on H/C/N/O, DZVP-MOLOPT-SR-GTH on Li, and cFIT3
auxiliary bases from BASIS_ADMM. The GPW cutoffs are 600/80 Ry; the inner SCF
tolerance is 10⁻⁷ Ha, and the Li population residual must be within 0.001 e.
ADMM uses the CP2K 2023.2 input keywords `METHOD BASIS_PROJECTION` and
`ADMM_PURIFICATION_METHOD MO_DIAG`, with `EXCH_CORRECTION_FUNC PBEX`.
The HFX interaction is truncated at 5 Å and checked against each cell size.
The orbital bases are atom-centered MOLOPT bases; a basis extension can be
assessed separately if the first energy comparison warrants it.

## Energies and interpretation

| JSON/CSV field | Definition | Interpretation |
|---|---|---|
| vertical_transfer_ev | Eᵦ(Rₐ) − Eₐ(Rₐ) | Transfer cost before structural relaxation |
| product_relaxation_stabilization_ev | Eᵦ(Rₐ) − Eᵦ(Rᵦ) | Energy gained by relaxing the electron-bearing state |
| relaxed_transfer_ev | Eᵦ(Rᵦ) − Eₐ(Rₐ) | Transfer cost after structural relaxation |
| reactant_reorganization_ev | Eₐ(Rᵦ) − Eₐ(Rₐ) | Reference-state reorganization energy |
| vertical_gap_at_product_geometry_ev | Eᵦ(Rᵦ) − Eₐ(Rᵦ) | State gap at the relaxed product geometry |

The main comparison is `eda_minus_tmeda_ev.relaxed_transfer_ev`. A negative
value means this modeled Li-to-solvent transfer has a lower relaxed energy cost
in EDA/THF. Compare energy differences within each composition before comparing
the amines. The output reports local-minimum electronic energies and the spread
between snapshots. This calculation covers Li electron transfer and solvent
relaxation; it contains no organic substrate or reaction pathway.

CP2K includes the residual cDFT term in its FORCE_EVAL energy. The analyzer
subtracts `(current_population − target_population) × constraint_strength`
from the final total, retaining both energies and the correction in the output.
The four single points must use matched coordinates, cells, and methods.
The same-state single points must also agree with the corresponding optimized
state energies within 0.02 eV.
If the product geometry has a lower reference-state energy than the optimized
reference local minimum, the record retains that observation as a diagnostic.
Unconverged optimizations or SCF calculations stop the workflow. A sign of the
EDA/TMEDA difference is never imposed as an acceptance condition.

## Existing source data

| System | Source campaign | Replica |
|---|---|---:|
| eda_3m | eda3m_pilot | 1 |
| tmeda_3m | amine_series_count_refinement | 1 |

For each source, the workflow reads the existing `analysis.json`,
`timeseries.csv`, composition spec, production TPR, and production XTC.
The TPR/XTC/spec hashes must match the completed classical analysis. Source
files are read-only inputs; these targets contain no MD rules. Adjust
`li_relaxation.source_campaigns` and `source_replicas` in `configs/methods.yaml`
if the completed files are under different campaigns.

Three snapshots per amine are selected from the equilibrated portion of the
existing trajectory using density and volume, separated by at least 1000 ps
or five estimated autocorrelation times. The bank records the chosen frames,
Li placement, cell volume, actual amine concentration, and one-Li concentration.
Snapshot spread comes from one source trajectory per amine.

## Run on TMC through Slurm

After updating the repository, start the pilot:

```bash
./run.sh dry-run --campaign li_eda_tmeda --target li_relaxation_pilot
./run.sh submit --campaign li_eda_tmeda --target li_relaxation_pilot
./run.sh inspect JOBID
```

`li_relaxation_pilot` uses one snapshot per amine: four optimizations and eight
single points. Its final summary contains the EDA/TMEDA energy comparison.
Each CP2K job uses eight MPI ranks and one OpenMP thread. The profiles retain
one active CP2K job and the twelve-CPU child-job budget, with a four-CPU controller.

The CP2K optimization wall-time request is 48 hours and single points request
24 hours. The controller's existing default is 12 hours. Set
`SOLVELEC_CONTROLLER_TIME` to a site-permitted duration sufficient for the
whole serial CP2K sequence before submitting longer runs, for example:

```bash
SOLVELEC_CONTROLLER_TIME=7-00:00:00 ./run.sh submit --campaign li_eda_tmeda --target li_relaxation_pilot
```

Once the pilot energies are available, the full three-snapshot comparison uses:

```bash
./run.sh submit --campaign li_eda_tmeda --target li_relaxation
```

The full target reuses completed pilot outputs and adds the remaining
snapshots. It contains twelve optimizations and twenty-four single points in
total. If a PBE/PBE0 energy comparison is useful, the optional benchmark adds
eight PBE single points on the pilot's PBE0 geometries, with no new optimization:

```bash
./run.sh submit --campaign li_eda_tmeda --target li_relaxation_benchmark
```

`li_relaxation_inputs` prepares the snapshot bank and the first-state inputs
without launching CP2K. `resume` can rerun an incomplete workflow using the
same target; a failed engine job starts a fresh attempt, retaining its previous
logs and CP2K restart files for inspection.

## Outputs

The default output root on the server is:

```text
/data/home/storage/Backup_Data/wangcx/li-thf-amine-solvated-electron/runs/li_eda_tmeda/li_relaxation
```

`pilot/summary.json` and `pilot/summary.csv` contain the initial comparison;
`full/summary.*` contains the three-snapshot means and sample standard
deviations. A single-snapshot standard deviation is `null`.

For each `SYSTEM/SNAPSHOT`, `opt/a/optimized.xyz` and
`opt/b/optimized.xyz` are the real optimized structures. Under
`sp/pbe0_admm/{aa,ba,bb,ab}/`, each job saves `electron.cube`, `spin.cube`,
`cp2k.out`, and `validation.json`. Total electron density includes all
valence electrons; spin density describes the alpha–beta difference.
Their coordinates and density grids can be used directly for subsequent
structure and density figures. Every attempt's raw files remain under its
`attempts/` directory.

## Implementation references

- [CP2K cDFT input documentation](https://manual.cp2k.org/cp2k-2023_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/QS/CDFT.html)
- [CP2K geometry optimization documentation](https://manual.cp2k.org/cp2k-2023_2-branch/CP2K_INPUT/MOTION/GEO_OPT.html)
- [CP2K 2023.2 cDFT energy term](https://github.com/cp2k/cp2k/blob/support/v2023.2/src/qs_cdft_methods.F)
- [CP2K 2023.2 total energy construction](https://github.com/cp2k/cp2k/blob/support/v2023.2/src/qs_ks_methods.F)
