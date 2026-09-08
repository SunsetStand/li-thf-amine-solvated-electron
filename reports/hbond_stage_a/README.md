# hbond_pilot Stage A report

This directory contains the shareable, Chinese-language report for the completed
`hbond_pilot` Stage A campaign:

- neat ethylenediamine (`pure_eda`, 64 EDA molecules);
- a 1:1 molar THF:EDA mixture (`eda_1to1`, 32 THF + 32 EDA molecules);
- three independent 20 ns production replicas per system.

The report is deliberately solvent-only. It contains no Li atom and no excess
electron. Its cavity is a geometric free-volume proxy, not an electron-density
or electron-localization result.

Open the finished report: [`hbond_stage_a_report_zh.pdf`](hbond_stage_a_report_zh.pdf).

## Rebuild

Install the report dependencies and run:

```bash
python -m pip install -e ".[report]"
python reports/hbond_stage_a/build_report.py
```

The builder reads only the compact, committed JSON/CSV/XYZ/PDB/SVG audit bundle
under `data/`; no trajectory, GROMACS installation, or scheduler is required.
It validates all six records, readiness gates, molecule counts, solvent-only
scope, frame counts, and SHA-256 provenance before writing any result.

Expected outputs:

- `hbond_stage_a_report_zh.pdf` - the illustrated report;
- `hbond_stage_a_metrics.json` - machine-readable aggregated metrics;
- `figures/*.png` - shareable plots and structure panels;
- `data/analysis/.../snapshot/cavity_local.pdb` and
  `cavity_hbonds.pml` - interactive 3D inspection inputs.

## Source calculation

The source campaign completed as Slurm controller job `17055` with 116/116
workflow steps. Publication metadata and file hashes are recorded in
`data/report_provenance.json`.
