# Stage B2C method-sensitivity benchmark

This target tests whether the small EDA-rich versus EDA-poor energy ordering from the
accepted three-replica Stage B2C production screen survives a tighter numerical setup and
a change from semilocal PBE to hybrid PBE0. It is a diagnostic benchmark, not a new
ensemble calculation and not a thermodynamic free energy.

## Immutable source and representative selection

The target requires the accepted files
`stage_b2c_preferential_production.summary.json` and
`stage_b2c_preferential_production.done`. It never recomputes their PBE baseline.
For each of `eda_1p5m` and `eda_3m`, the plan step selects and verifies two distinct
replicas:

- the largest positive EDA-rich minus EDA-poor attachment-energy proxy;
- the smallest (most negative) EDA-rich minus EDA-poor proxy.

The checked configuration currently corresponds to EDA 1.5 M replicas 2 and 3 and EDA
3 M replicas 1 and 3. If the accepted production summary no longer yields these extrema,
planning fails instead of silently benchmarking different structures.

## Electronic-structure matrix

Every selected replica retains both existing seed roles (`eda_rich`, `eda_poor`) and the
matched neutral/anion fixed-nuclei pair. Two new methods are evaluated:

- `pbe_converged`: PBE-D3(BJ), TZV2P-MOLOPT-GTH, 600/80 Ry, EPS_SCF 1e-7;
- `pbe0_admm`: PBE0-D3(BJ), the same orbital basis and grid, cFIT3 ADMM, 25% exact
  exchange, and a 5 angstrom truncated-exchange cutoff.

This gives 4 selected system/replica pairs x 2 seed roles x 2 charge states x 2 methods,
or 32 new CP2K single points. The complete DAG contains 67 jobs: one plan, 16 renders,
32 CP2K runs, 16 analyses, one summary, and one gate. It contains no GROMACS jobs.

## Run on TMC

After updating the server checkout, preview once:

```bash
./run.sh dry-run --campaign pilot --target stage_b2c_method_benchmark
```

The expected dry-run total is 67 and the rule table must show 32 CP2K runs. Submit with:

```bash
./run.sh submit --campaign pilot --target stage_b2c_method_benchmark
```

The number 32 above is the number of independent single-point calculations,
not the CPU allocation per calculation. On TMC each CP2K job requests twelve
single-threaded MPI ranks. The `cp2k_slots=1` and `cpu_slots=12` guards prevent
other child jobs from overlapping it; with the four-CPU controller, total
workflow allocation is at most sixteen CPUs. Do not override these guards or
raise the rank count without a representative scaling benchmark and site
approval.

Inspect the controller using `./run.sh inspect JOBID`. Success requires both:

- `runs/pilot/stage_b2c_method_benchmark.summary.json` with `ready: true`;
- `runs/pilot/stage_b2c_method_benchmark.done` whose SHA-256 matches the summary.

On TMC these paths are rooted below the configured storage run directory.

## Interpretation

The summary reports, for every selected pair and method, the raw EDA-rich minus EDA-poor
energy difference, the shift from the accepted PBE baseline, the 0.05 eV preference class,
and compact spin-localization diagnostics. A method-dependent sign, class, or shift is a
scientific result (`method_sensitive_or_unresolved`), not a workflow failure.

Even agreement cannot establish preferential-solvation thermodynamics: the test uses only
selected extrema, fixed nuclei, uncorrected charged periodic cells, and PBE0/ADMM rather
than a complete-basis reference. It addresses whether the current ordering is robust enough
to justify more sampling or higher-level finite-size work.
