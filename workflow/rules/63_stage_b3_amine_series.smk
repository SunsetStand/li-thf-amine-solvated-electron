# Prepare only original amine sources still selected by Stage B3. All raw
# trajectory paths are immutable params: this target has no MD producer edges.
STAGE_B3_ORIGINAL_PAIRS = [
    pair for pair in STAGE_B3_ENVIRONMENT_PAIRS
    if stage_b3_source_campaign(pair[0]) == "amine_series_pilot"
]
STAGE_B3_SOURCE_HANDOFF_ROOT = f"{RUN_ROOT}/pilot/stage_b3/source_handoff/amine_series_pilot"


rule prepare_stage_b3_source_handoff:
    input:
        methods="configs/methods.yaml",
        scripts=[STAGE_B3_AMINE_SERIES, ANALYZE_CLASSICAL, PREPARE_STAGE_B],
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        manifest=f"{STAGE_B3_SOURCE_HANDOFF_ROOT}/{{system}}/r{{replica}}/candidates/manifest.json",
        summary=f"{STAGE_B3_SOURCE_HANDOFF_ROOT}/{{system}}/r{{replica}}/stage_b_candidates.summary.json",
        gate=f"{STAGE_B3_SOURCE_HANDOFF_ROOT}/{{system}}/r{{replica}}/stage_b_candidates.done",
    params:
        source_campaign=lambda wc: stage_b3_source_campaign(wc.system),
        spec=lambda wc: stage_b3_artifact(f"{RUN_ROOT}/amine_series_pilot/specs/{wc.system}/r{wc.replica}.json"),
        validation=lambda wc: stage_b3_artifact(f"{RUN_ROOT}/amine_series_pilot/classical/{wc.system}/r{wc.replica}/pilot/validation.json"),
        tpr=lambda wc: stage_b3_artifact(f"{RUN_ROOT}/amine_series_pilot/classical/{wc.system}/r{wc.replica}/pilot/production/production.tpr"),
        trajectory=lambda wc: stage_b3_artifact(f"{RUN_ROOT}/amine_series_pilot/classical/{wc.system}/r{wc.replica}/pilot/production/production.xtc"),
        directory=lambda wc: f"{STAGE_B3_SOURCE_HANDOFF_ROOT}/{wc.system}/r{wc.replica}",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=180,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {STAGE_B3_AMINE_SERIES:q} source-handoff "
        "--campaign {CAMPAIGN:q} --source-campaign {params.source_campaign:q} "
        "--spec {params.spec:q} --validation {params.validation:q} --tpr {params.tpr:q} "
        "--trajectory {params.trajectory:q} --methods {input.methods:q} --output-dir {params.directory:q}"


rule stage_b3_source_handoff:
    input:
        [f"{STAGE_B3_SOURCE_HANDOFF_ROOT}/{system}/r{replica}/stage_b_candidates.done"
         for system, replica in STAGE_B3_ORIGINAL_PAIRS]
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b3_source_handoff.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "touch {output:q}"


rule analyze_stage_b3_existing_nh:
    input:
        methods="configs/methods.yaml",
        systems="configs/systems.yaml",
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/reanalysis/"
            "{system}/r{replica}/{seed_role}.json"
        )
    params:
        production_gate=stage_b3_b2c_production_gate,
        analysis=lambda wildcards: stage_b3_existing_analysis(
            wildcards.system, wildcards.replica, wildcards.seed_role
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed_role="eda_(rich|poor)",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} reanalyze "
        "--analysis {params.analysis:q} --methods {input.methods:q} "
        "--systems {input.systems:q} --output {output:q}"


rule summarize_stage_b3_nh_reanalysis:
    input:
        records=STAGE_B3_REANALYSIS_RECORDS,
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b3_nh_reanalysis.summary.json"
    params:
        expected_count=len(STAGE_B3_REANALYSIS_RECORDS),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} "
        "reanalysis-summary --campaign {CAMPAIGN:q} "
        "--expected-count {params.expected_count} --records {input.records:q} "
        "--output {output:q}"


rule stage_b3_nh_reanalysis:
    input:
        summary=rules.summarize_stage_b3_nh_reanalysis.output,
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b3_nh_reanalysis.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"


rule screen_stage_b3_amine_environment:
    input:
        methods="configs/methods.yaml",
        systems="configs/systems.yaml",
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/environment/{{system}}/r{{replica}}.json"
    params:
        source_campaign=lambda wildcards: stage_b3_source_campaign(wildcards.system),
        candidate_gate=stage_b3_candidate_gate,
        candidate_summary=lambda wildcards: stage_b3_artifact(
            stage_b3_candidate_path(wildcards, "stage_b_candidates.summary.json")
        ),
        candidate_manifest=lambda wildcards: stage_b3_artifact(
            stage_b3_candidate_path(wildcards, "manifest")
        ),
        spec=lambda wildcards: stage_b3_artifact(
            f"{RUN_ROOT}/{stage_b3_source_campaign(wildcards.system)}/specs/"
            f"{wildcards.system}/r{wildcards.replica}.json"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} screen "
        "--campaign {CAMPAIGN:q} --source-campaign {params.source_campaign:q} "
        "--candidate-manifest {params.candidate_manifest:q} "
        "--candidate-summary {params.candidate_summary:q} "
        "--candidate-gate {params.candidate_gate:q} --spec {params.spec:q} "
        "--methods {input.methods:q} --systems {input.systems:q} --output {output:q}"


rule summarize_stage_b3_amine_environment:
    input:
        records=STAGE_B3_ENVIRONMENT_RECORDS,
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b3_amine_environment.summary.json"
    params:
        systems=" ".join(shlex.quote(value) for value in STAGE_B3_ENVIRONMENT_SYSTEMS),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} "
        "environment-summary --campaign {CAMPAIGN:q} "
        "--expected-systems {params.systems} --records {input.records:q} "
        "--output {output:q}"


rule stage_b3_amine_environment:
    input:
        summary=rules.summarize_stage_b3_amine_environment.output,
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b3_amine_environment.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"


rule prepare_stage_b3_electronic_seeds:
    input:
        screen=rules.screen_stage_b3_amine_environment.output,
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        manifest=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/"
            "{system}/r{replica}/manifest.json"
        ),
        coordinates=[
            (
                f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/"
                f"{{system}}/r{{replica}}/{seed_role}/coordinates.xyz"
            )
            for seed_role in STAGE_B3_SEED_ROLES
        ],
        cells=[
            (
                f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/"
                f"{{system}}/r{{replica}}/{seed_role}/cell.inc"
            )
            for seed_role in STAGE_B3_SEED_ROLES
        ],
        metadata=[
            (
                f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/"
                f"{{system}}/r{{replica}}/{seed_role}/metadata.json"
            )
            for seed_role in STAGE_B3_SEED_ROLES
        ],
    params:
        candidate_manifest=lambda wildcards: stage_b3_artifact(
            stage_b3_candidate_path(wildcards, "manifest")
        ),
        output_dir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/"
            f"{wildcards.system}/r{wildcards.replica}"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} prepare "
        "--screen {input.screen:q} --candidate-manifest {params.candidate_manifest:q} "
        "--output-dir {params.output_dir:q} --output {output.manifest:q}"


rule render_stage_b3_electronic_pair:
    input:
        manifest=rules.prepare_stage_b3_electronic_seeds.output.manifest,
        methods="configs/methods.yaml",
        template="workflow/templates/cp2k/stage_b2c_preferential_smoke.inp.tpl",
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        neutral=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/"
            "{system}/r{replica}/{seed_role}/neutral/cp2k.inp"
        ),
        anion=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/"
            "{system}/r{replica}/{seed_role}/anion/cp2k.inp"
        ),
    params:
        project_prefix=lambda wildcards: (
            f"{wildcards.system}_r{wildcards.replica}_{wildcards.seed_role}"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed_role="nh_(facing|control)",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} render "
        "--manifest {input.manifest:q} --methods {input.methods:q} "
        "--template {input.template:q} --seed-role {wildcards.seed_role:q} "
        "--project-prefix {params.project_prefix:q} --neutral-output {output.neutral:q} "
        "--anion-output {output.anion:q}"


rule run_stage_b3_electronic_neutral:
    input:
        cp2k=rules.render_stage_b3_electronic_pair.output.neutral,
        engine_runner=ENGINE_RUNNER,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/"
            "{system}/r{replica}/{seed_role}/neutral/cp2k.out"
        )
    params:
        workdir=lambda wildcards: stage_b3_directory(
            wildcards.system, wildcards.replica, wildcards.seed_role, "neutral"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed_role="nh_(facing|control)",
    threads: 1
    resources:
        tasks=8,
        cp2k_slots=1,
        cpu_slots=8,
        mpi="mpirun",
        mem_mb=128000,
        runtime=720,
    shell:
        "bash {STAGE_RUNNER:q} cdft -- {PYTHON} {input.engine_runner:q} "
        "--engine cp2k --output {output:q} --cwd {params.workdir:q} -- "
        "{resources.mpi} -n {resources.tasks} cp2k.psmp -i {input.cp2k:q}"


rule run_stage_b3_electronic_anion:
    input:
        cp2k=rules.render_stage_b3_electronic_pair.output.anion,
        engine_runner=ENGINE_RUNNER,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        cp2k=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/"
            "{system}/r{replica}/{seed_role}/anion/cp2k.out"
        ),
        electron_cube=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/{{system}}/r{{replica}}/"
            "{seed_role}/anion/{system}_r{replica}_{seed_role}_anion_"
            "stage_b3-ELECTRON_DENSITY-1_0.cube"
        ),
        spin_cube=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/{{system}}/r{{replica}}/"
            "{seed_role}/anion/{system}_r{replica}_{seed_role}_anion_"
            "stage_b3-SPIN_DENSITY-1_0.cube"
        ),
    params:
        workdir=lambda wildcards: stage_b3_directory(
            wildcards.system, wildcards.replica, wildcards.seed_role, "anion"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed_role="nh_(facing|control)",
    threads: 1
    resources:
        tasks=8,
        cp2k_slots=1,
        cpu_slots=8,
        mpi="mpirun",
        mem_mb=128000,
        runtime=720,
    shell:
        "bash {STAGE_RUNNER:q} cdft -- {PYTHON} {input.engine_runner:q} "
        "--engine cp2k --output {output.cp2k:q} --cwd {params.workdir:q} -- "
        "{resources.mpi} -n {resources.tasks} cp2k.psmp -i {input.cp2k:q}"


rule analyze_stage_b3_electronic_pair:
    input:
        manifest=rules.prepare_stage_b3_electronic_seeds.output.manifest,
        neutral_input=rules.render_stage_b3_electronic_pair.output.neutral,
        anion_input=rules.render_stage_b3_electronic_pair.output.anion,
        neutral_output=rules.run_stage_b3_electronic_neutral.output,
        anion_output=rules.run_stage_b3_electronic_anion.output.cp2k,
        electron_cube=rules.run_stage_b3_electronic_anion.output.electron_cube,
        spin_cube=rules.run_stage_b3_electronic_anion.output.spin_cube,
        methods="configs/methods.yaml",
        systems="configs/systems.yaml",
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b3/electronic/"
            "{system}/r{replica}/{seed_role}/analysis.json"
        )
    params:
        spec=lambda wildcards: stage_b3_artifact(
            f"{RUN_ROOT}/{stage_b3_source_campaign(wildcards.system)}/specs/"
            f"{wildcards.system}/r{wildcards.replica}.json"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed_role="nh_(facing|control)",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} analyze "
        "--manifest {input.manifest:q} --seed-role {wildcards.seed_role:q} "
        "--spec {params.spec:q} --methods {input.methods:q} --systems {input.systems:q} "
        "--neutral-input {input.neutral_input:q} --neutral-output {input.neutral_output:q} "
        "--anion-input {input.anion_input:q} --anion-output {input.anion_output:q} "
        "--spin-cube {input.spin_cube:q} --electron-cube {input.electron_cube:q} "
        "--output {output:q}"


rule summarize_stage_b3_amine_pilot:
    input:
        records=STAGE_B3_PILOT_ANALYSES,
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b3_amine_pilot.summary.json"
    params:
        systems=" ".join(shlex.quote(value) for value in STAGE_B3_PILOT_SYSTEMS),
        replicas=" ".join(str(value) for value in sorted({r for _s, r in STAGE_B3_PILOT_PAIRS})),
        tolerance=float(STAGE_B3_SETTINGS["energy_tie_tolerance_ev"]),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} "
        "pilot-summary --campaign {CAMPAIGN:q} --expected-systems {params.systems} "
        "--expected-replicas {params.replicas} --energy-tolerance-ev {params.tolerance} "
        "--records {input.records:q} --output {output:q}"


rule stage_b3_amine_pilot:
    input:
        summary=rules.summarize_stage_b3_amine_pilot.output,
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b3_amine_pilot.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"


rule summarize_stage_b3_amine_series:
    input:
        reanalysis_summary=rules.summarize_stage_b3_nh_reanalysis.output,
        reanalysis_gate=rules.stage_b3_nh_reanalysis.output,
        environment_summary=rules.summarize_stage_b3_amine_environment.output,
        environment_gate=rules.stage_b3_amine_environment.output,
        pilot_summary=rules.summarize_stage_b3_amine_pilot.output,
        pilot_gate=rules.stage_b3_amine_pilot.output,
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b3_amine_series.summary.json"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} combine "
        "--campaign {CAMPAIGN:q} --reanalysis-summary {input.reanalysis_summary:q} "
        "--reanalysis-gate {input.reanalysis_gate:q} "
        "--environment-summary {input.environment_summary:q} "
        "--environment-gate {input.environment_gate:q} "
        "--pilot-summary {input.pilot_summary:q} --pilot-gate {input.pilot_gate:q} "
        "--output {output:q}"


rule stage_b3_amine_series:
    input:
        summary=rules.summarize_stage_b3_amine_series.output,
        script=STAGE_B3_AMINE_SERIES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b3_amine_series.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"
