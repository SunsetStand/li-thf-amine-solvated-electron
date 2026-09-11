# Complete Stage-B report

`build_report.py` assembles the accepted Stage-B pilot artifacts into a Chinese
PDF, five figures, compact audit data, machine-readable metrics, and a full
SHA-256 provenance manifest. It reads existing candidate, CP2K, and cube files
in place; it never runs or modifies an electronic-structure or MD calculation.

On TMC-AMD, use the report-only Snakemake target after
`stage_b_localization.done` has passed:

```bash
./run.sh dry-run --campaign pilot --target stage_b_report
./run.sh submit --campaign pilot --target stage_b_report
```

The dry-run should contain exactly two lightweight rules:

- `build_stage_b_report` (one report assembly job);
- `stage_b_report` (one SHA-256 completion gate).

It must not contain CP2K or GROMACS rules. Outputs are written to
`$SOLVELEC_RUN_ROOT/pilot/stage_b_report/`; the final gate is
`$SOLVELEC_RUN_ROOT/pilot/stage_b_report.done`.

The report deliberately labels the current result as a completed numerical
pilot rather than a production mechanism or stability result. Large cube and
CP2K output files remain in storage and are referenced by path, size, and hash;
only compact summaries, candidate geometry, CP2K inputs, figures, and the PDF
are copied into the report directory.
