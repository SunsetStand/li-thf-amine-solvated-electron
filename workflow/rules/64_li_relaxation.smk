from solvelec.li_relaxation import settings_from as li_settings_from

LI_SETTINGS = li_settings_from(METHOD_CATALOG)
LI_SCRIPT = str(ROOT / "workflow/scripts/li_relaxation.py")
LI_BASE = f"{RUN_ROOT}/{CAMPAIGN}/li_relaxation"
LI_SNAPSHOTS = [f"s{i:02d}" for i in range(1, int(LI_SETTINGS["snapshots_per_system"])+1)]
LI_SYSTEMS = LI_SETTINGS["systems"]
LI_VARIANT = LI_SETTINGS["optimization_variant"]
LI_SCOPES = {
    "pilot": (LI_SNAPSHOTS[:1], [LI_VARIANT]),
    "full": (LI_SNAPSHOTS, [LI_VARIANT]),
    "benchmark": (LI_SNAPSHOTS[:1], list(LI_SETTINGS["variants"])),
}


def li_require_campaign():
    if CAMPAIGN != LI_SETTINGS["campaign"]:
        raise ValueError(f"Li relaxation requires --campaign {LI_SETTINGS['campaign']}")


def li_source(wildcards, kind):
    li_require_campaign()
    source = LI_SETTINGS["source_campaigns"][wildcards.system]
    replica = int(LI_SETTINGS["source_replicas"][wildcards.system])
    root = Path(RUN_ROOT) / source
    paths = {
        "analysis": root / f"analysis/{wildcards.system}/r{replica}/analysis.json",
        "timeseries": root / f"analysis/{wildcards.system}/r{replica}/timeseries.csv",
        "spec": root / f"specs/{wildcards.system}/r{replica}.json",
        "tpr": root / f"classical/{wildcards.system}/r{replica}/pilot/production/production.tpr",
        "trajectory": root / f"classical/{wildcards.system}/r{replica}/pilot/production/production.xtc",
    }
    if not paths[kind].is_file():
        raise ValueError(f"Existing MD handoff missing: {paths[kind]}; this target will not rerun MD")
    return str(paths[kind])


def li_opt_coordinates(wildcards):
    if wildcards.state == "a":
        return f"{LI_BASE}/{wildcards.system}/seeds/{wildcards.snapshot}/coordinates.xyz"
    return f"{LI_BASE}/{wildcards.system}/{wildcards.snapshot}/opt/a/optimized.xyz"


def li_opt_validation(wildcards):
    if wildcards.state == "a":
        return []
    return [f"{LI_BASE}/{wildcards.system}/{wildcards.snapshot}/opt/a/validation.json"]


def li_matrix_records(wildcards):
    snapshots, variants = LI_SCOPES[wildcards.scope]
    return [f"{LI_BASE}/{system}/{snapshot}/analysis/{variant}.json"
            for system in LI_SYSTEMS for snapshot in snapshots for variant in variants]


rule prepare_li_relaxation_seeds:
    input:
        analysis=lambda w: ancient(li_source(w, "analysis")),
        timeseries=lambda w: ancient(li_source(w, "timeseries")),
        spec=lambda w: ancient(li_source(w, "spec")),
        tpr=lambda w: ancient(li_source(w, "tpr")),
        trajectory=lambda w: ancient(li_source(w, "trajectory")),
        methods="configs/methods.yaml",
        script=LI_SCRIPT,
        sources=SOLVELEC_SOURCES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        bank=f"{LI_BASE}/{{system}}/seeds/bank.json",
        xyz=[f"{LI_BASE}/{{system}}/seeds/{snapshot}/coordinates.xyz" for snapshot in LI_SNAPSHOTS],
        cells=[f"{LI_BASE}/{{system}}/seeds/{snapshot}/cell.inc" for snapshot in LI_SNAPSHOTS],
    wildcard_constraints:
        system="eda_3m|tmeda_3m",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        directory=lambda w: f"{LI_BASE}/{w.system}/seeds",
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} prepare "
        "--run-root {RUN_ROOT:q} --methods {input.methods:q} --system {wildcards.system:q} "
        "--output-dir {params.directory:q}"


rule render_li_relaxation_opt:
    input:
        bank=f"{LI_BASE}/{{system}}/seeds/bank.json",
        coordinates=li_opt_coordinates,
        validation=li_opt_validation,
        methods="configs/methods.yaml",
        template="workflow/templates/cp2k/li_relaxation.inp.tpl",
        script=LI_SCRIPT,
        sources=SOLVELEC_SOURCES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        cp2k=f"{LI_BASE}/{{system}}/{{snapshot}}/opt/{{state}}/cp2k.inp",
        job=f"{LI_BASE}/{{system}}/{{snapshot}}/opt/{{state}}/cp2k.job.json",
    wildcard_constraints:
        system="eda_3m|tmeda_3m", snapshot="|".join(LI_SNAPSHOTS), state="a|b",
    params:
        variant=LI_VARIANT,
        geometry=lambda w: "seed" if w.state == "a" else "a",
        validation=lambda w, input: (
            "--geometry-validation " + shlex.quote(str(input.validation[0]))
            if input.validation else ""
        ),
    threads: 4
    resources:
        mem_mb=4000, runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} render "
        "--methods {input.methods:q} --bank {input.bank:q} --snapshot {wildcards.snapshot:q} "
        "--state {wildcards.state:q} --variant {params.variant:q} --optimize "
        "--geometry-state {params.geometry:q} {params.validation} "
        "--coordinates {input.coordinates:q} --template {input.template:q} --output {output.cp2k:q}"


rule run_li_relaxation_opt:
    input:
        cp2k=rules.render_li_relaxation_opt.output.cp2k,
        job=rules.render_li_relaxation_opt.output.job,
        script=LI_SCRIPT,
        sources=SOLVELEC_SOURCES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        cp2k=f"{LI_BASE}/{{system}}/{{snapshot}}/opt/{{state}}/cp2k.out",
        xyz=f"{LI_BASE}/{{system}}/{{snapshot}}/opt/{{state}}/optimized.xyz",
        validation=f"{LI_BASE}/{{system}}/{{snapshot}}/opt/{{state}}/validation.json",
    wildcard_constraints:
        system="eda_3m|tmeda_3m", snapshot="|".join(LI_SNAPSHOTS), state="a|b",
    params:
        directory=lambda w: f"{LI_BASE}/{w.system}/{w.snapshot}/opt/{w.state}",
    threads: 1
    resources:
        tasks=8, cp2k_slots=1, cpu_slots=8, mpi="mpirun", mem_mb=128000, runtime=2880,
    shell:
        "bash {STAGE_RUNNER:q} cdft -- {PYTHON} {input.script:q} run "
        "--job {input.job:q} --workdir {params.directory:q} -- "
        "{resources.mpi} -n {resources.tasks} cp2k.psmp"


rule render_li_relaxation_sp:
    input:
        bank=f"{LI_BASE}/{{system}}/seeds/bank.json",
        coordinates=lambda w: f"{LI_BASE}/{w.system}/{w.snapshot}/opt/{w.pair[1]}/optimized.xyz",
        validation=lambda w: f"{LI_BASE}/{w.system}/{w.snapshot}/opt/{w.pair[1]}/validation.json",
        methods="configs/methods.yaml",
        template="workflow/templates/cp2k/li_relaxation.inp.tpl",
        script=LI_SCRIPT,
        sources=SOLVELEC_SOURCES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        cp2k=f"{LI_BASE}/{{system}}/{{snapshot}}/sp/{{variant}}/{{pair}}/cp2k.inp",
        job=f"{LI_BASE}/{{system}}/{{snapshot}}/sp/{{variant}}/{{pair}}/cp2k.job.json",
    wildcard_constraints:
        system="eda_3m|tmeda_3m", snapshot="|".join(LI_SNAPSHOTS),
        variant="|".join(LI_SETTINGS["variants"]), pair="aa|ba|bb|ab",
    params:
        state=lambda w: w.pair[0], geometry=lambda w: w.pair[1],
    threads: 4
    resources:
        mem_mb=4000, runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} render "
        "--methods {input.methods:q} --bank {input.bank:q} --snapshot {wildcards.snapshot:q} "
        "--state {params.state:q} --variant {wildcards.variant:q} "
        "--geometry-state {params.geometry:q} --geometry-validation {input.validation:q} "
        "--coordinates {input.coordinates:q} --template {input.template:q} --output {output.cp2k:q}"


rule run_li_relaxation_sp:
    input:
        cp2k=rules.render_li_relaxation_sp.output.cp2k,
        job=rules.render_li_relaxation_sp.output.job,
        script=LI_SCRIPT,
        sources=SOLVELEC_SOURCES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        cp2k=f"{LI_BASE}/{{system}}/{{snapshot}}/sp/{{variant}}/{{pair}}/cp2k.out",
        validation=f"{LI_BASE}/{{system}}/{{snapshot}}/sp/{{variant}}/{{pair}}/validation.json",
        electron=f"{LI_BASE}/{{system}}/{{snapshot}}/sp/{{variant}}/{{pair}}/electron.cube",
        spin=f"{LI_BASE}/{{system}}/{{snapshot}}/sp/{{variant}}/{{pair}}/spin.cube",
    wildcard_constraints:
        system="eda_3m|tmeda_3m", snapshot="|".join(LI_SNAPSHOTS),
        variant="|".join(LI_SETTINGS["variants"]), pair="aa|ba|bb|ab",
    params:
        directory=lambda w: f"{LI_BASE}/{w.system}/{w.snapshot}/sp/{w.variant}/{w.pair}",
    threads: 1
    resources:
        tasks=8, cp2k_slots=1, cpu_slots=8, mpi="mpirun", mem_mb=128000, runtime=1440,
    shell:
        "bash {STAGE_RUNNER:q} cdft -- {PYTHON} {input.script:q} run "
        "--job {input.job:q} --workdir {params.directory:q} -- "
        "{resources.mpi} -n {resources.tasks} cp2k.psmp"


rule analyze_li_relaxation_matrix:
    input:
        validations=lambda w: [f"{LI_BASE}/{w.system}/{w.snapshot}/sp/{w.variant}/{pair}/validation.json"
                               for pair in ("aa", "ba", "bb", "ab")],
        bank=f"{LI_BASE}/{{system}}/seeds/bank.json",
        methods="configs/methods.yaml",
        script=LI_SCRIPT,
        sources=SOLVELEC_SOURCES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{LI_BASE}/{{system}}/{{snapshot}}/analysis/{{variant}}.json",
    wildcard_constraints:
        system="eda_3m|tmeda_3m", snapshot="|".join(LI_SNAPSHOTS),
        variant="|".join(LI_SETTINGS["variants"]),
    params:
        directory=lambda w: f"{LI_BASE}/{w.system}/{w.snapshot}",
    threads: 4
    resources:
        mem_mb=4000, runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} analyze "
        "--methods {input.methods:q} --bank {input.bank:q} --snapshot {wildcards.snapshot:q} "
        "--variant {wildcards.variant:q} --snapshot-dir {params.directory:q} --output {output:q}"


rule summarize_li_relaxation:
    input:
        records=li_matrix_records,
        methods="configs/methods.yaml",
        script=LI_SCRIPT,
        sources=SOLVELEC_SOURCES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        summary=f"{LI_BASE}/{{scope}}/summary.json", csv=f"{LI_BASE}/{{scope}}/summary.csv",
    wildcard_constraints:
        scope="pilot|full|benchmark",
    params:
        snapshots=lambda w: " ".join(shlex.quote(x) for x in LI_SCOPES[w.scope][0]),
        variants=lambda w: " ".join(shlex.quote(x) for x in LI_SCOPES[w.scope][1]),
    threads: 4
    resources:
        mem_mb=4000, runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} summary "
        "--campaign {CAMPAIGN:q} --methods {input.methods:q} --snapshots {params.snapshots} "
        "--variants {params.variants} --records {input.records:q} --output {output.summary:q}"


rule validate_li_relaxation_summary:
    input:
        summary=f"{LI_BASE}/{{scope}}/summary.json",
        script=LI_SCRIPT,
        sources=SOLVELEC_SOURCES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{LI_BASE}/{{scope}}/done",
    wildcard_constraints:
        scope="pilot|full|benchmark",
    threads: 4
    resources:
        mem_mb=4000, runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"


rule li_relaxation_inputs:
    input:
        expand(f"{LI_BASE}/{{system}}/seeds/bank.json", system=LI_SYSTEMS),
        expand(f"{LI_BASE}/{{system}}/s01/opt/a/cp2k.inp", system=LI_SYSTEMS),


rule li_relaxation_pilot:
    input:
        f"{LI_BASE}/pilot/done",


rule li_relaxation_benchmark:
    input:
        f"{LI_BASE}/benchmark/done",


rule li_relaxation:
    input:
        f"{LI_BASE}/full/done",
