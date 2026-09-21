rule plan_stage_b2c_method_benchmark:
    input:
        methods="configs/methods.yaml",
        script=STAGE_B2C_METHOD_BENCHMARK,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_method_benchmark/plan.json"
    params:
        production_summary=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential_production.summary.json"
        ),
        production_gate=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential_production.done"
        ),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} plan "
        "--campaign {CAMPAIGN:q} --production-summary {params.production_summary:q} "
        "--production-gate {params.production_gate:q} --methods {input.methods:q} "
        "--output {output:q}"


rule render_stage_b2c_method_benchmark_pair:
    input:
        plan=rules.plan_stage_b2c_method_benchmark.output,
        methods="configs/methods.yaml",
        template="workflow/templates/cp2k/stage_b2c_method_benchmark.inp.tpl",
        script=STAGE_B2C_METHOD_BENCHMARK,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        neutral=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_method_benchmark/"
            "{system}/r{replica}/{seed_role}/{method_variant}/neutral/cp2k.inp"
        ),
        anion=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_method_benchmark/"
            "{system}/r{replica}/{seed_role}/{method_variant}/anion/cp2k.inp"
        ),
    params:
        manifest=lambda wildcards: stage_b2c_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/{wildcards.system}/"
            f"r{wildcards.replica}/manifest.json"
        ),
        project_prefix=lambda wildcards: (
            f"{wildcards.system}_r{wildcards.replica}_{wildcards.seed_role}_"
            f"{wildcards.method_variant}"
        ),
    wildcard_constraints:
        system="(eda_1p5m|eda_3m)",
        replica="[1-9][0-9]*",
        seed_role="eda_(rich|poor)",
        method_variant="(pbe_converged|pbe0_admm)",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} render "
        "--system {wildcards.system:q} --replica {wildcards.replica} "
        "--seed-role {wildcards.seed_role:q} --method-variant {wildcards.method_variant:q} "
        "--plan {input.plan:q} --manifest {params.manifest:q} --methods {input.methods:q} "
        "--template {input.template:q} --project-prefix {params.project_prefix:q} "
        "--neutral-output {output.neutral:q} --anion-output {output.anion:q}"


rule run_stage_b2c_method_benchmark_neutral:
    input:
        cp2k=rules.render_stage_b2c_method_benchmark_pair.output.neutral,
        engine_runner=ENGINE_RUNNER,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_method_benchmark/"
            "{system}/r{replica}/{seed_role}/{method_variant}/neutral/cp2k.out"
        )
    params:
        workdir=lambda wildcards: stage_b2c_method_benchmark_directory(
            wildcards.system,
            wildcards.replica,
            wildcards.seed_role,
            wildcards.method_variant,
            "neutral",
        ),
    wildcard_constraints:
        system="(eda_1p5m|eda_3m)",
        replica="[1-9][0-9]*",
        seed_role="eda_(rich|poor)",
        method_variant="(pbe_converged|pbe0_admm)",
    threads: 1
    resources:
        tasks=12,
        cp2k_slots=1,
        cpu_slots=12,
        mpi="mpirun",
        mem_mb=128000,
        runtime=1440,
    shell:
        "bash {STAGE_RUNNER:q} cdft -- {PYTHON} {input.engine_runner:q} "
        "--engine cp2k --output {output:q} --cwd {params.workdir:q} -- "
        "{resources.mpi} -n {resources.tasks} cp2k.psmp -i {input.cp2k:q}"


rule run_stage_b2c_method_benchmark_anion:
    input:
        cp2k=rules.render_stage_b2c_method_benchmark_pair.output.anion,
        engine_runner=ENGINE_RUNNER,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        cp2k=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_method_benchmark/"
            "{system}/r{replica}/{seed_role}/{method_variant}/anion/cp2k.out"
        ),
        electron_cube=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_method_benchmark/{{system}}/r{{replica}}/"
            "{seed_role}/{method_variant}/anion/{system}_r{replica}_{seed_role}_"
            "{method_variant}_anion_stage_b2c_method_benchmark-ELECTRON_DENSITY-1_0.cube"
        ),
        spin_cube=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_method_benchmark/{{system}}/r{{replica}}/"
            "{seed_role}/{method_variant}/anion/{system}_r{replica}_{seed_role}_"
            "{method_variant}_anion_stage_b2c_method_benchmark-SPIN_DENSITY-1_0.cube"
        ),
    params:
        workdir=lambda wildcards: stage_b2c_method_benchmark_directory(
            wildcards.system,
            wildcards.replica,
            wildcards.seed_role,
            wildcards.method_variant,
            "anion",
        ),
    wildcard_constraints:
        system="(eda_1p5m|eda_3m)",
        replica="[1-9][0-9]*",
        seed_role="eda_(rich|poor)",
        method_variant="(pbe_converged|pbe0_admm)",
    threads: 1
    resources:
        tasks=12,
        cp2k_slots=1,
        cpu_slots=12,
        mpi="mpirun",
        mem_mb=128000,
        runtime=1440,
    shell:
        "bash {STAGE_RUNNER:q} cdft -- {PYTHON} {input.engine_runner:q} "
        "--engine cp2k --output {output.cp2k:q} --cwd {params.workdir:q} -- "
        "{resources.mpi} -n {resources.tasks} cp2k.psmp -i {input.cp2k:q}"


rule analyze_stage_b2c_method_benchmark_pair:
    input:
        plan=rules.plan_stage_b2c_method_benchmark.output,
        neutral_input=rules.render_stage_b2c_method_benchmark_pair.output.neutral,
        anion_input=rules.render_stage_b2c_method_benchmark_pair.output.anion,
        neutral_output=rules.run_stage_b2c_method_benchmark_neutral.output,
        anion_output=rules.run_stage_b2c_method_benchmark_anion.output.cp2k,
        electron_cube=rules.run_stage_b2c_method_benchmark_anion.output.electron_cube,
        spin_cube=rules.run_stage_b2c_method_benchmark_anion.output.spin_cube,
        methods="configs/methods.yaml",
        systems="configs/systems.yaml",
        script=STAGE_B2C_METHOD_BENCHMARK,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_method_benchmark/"
            "{system}/r{replica}/{seed_role}/{method_variant}/analysis.json"
        )
    params:
        manifest=lambda wildcards: stage_b2c_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/{wildcards.system}/"
            f"r{wildcards.replica}/manifest.json"
        ),
        spec=lambda wildcards: stage_b2c_artifact(
            f"{RUN_ROOT}/{stage_b2c_source_campaign(wildcards.system)}/specs/"
            f"{wildcards.system}/r{wildcards.replica}.json"
        ),
    wildcard_constraints:
        system="(eda_1p5m|eda_3m)",
        replica="[1-9][0-9]*",
        seed_role="eda_(rich|poor)",
        method_variant="(pbe_converged|pbe0_admm)",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} analyze "
        "--campaign {CAMPAIGN:q} --system {wildcards.system:q} "
        "--replica {wildcards.replica} --seed-role {wildcards.seed_role:q} "
        "--method-variant {wildcards.method_variant:q} --plan {input.plan:q} "
        "--manifest {params.manifest:q} --spec {params.spec:q} "
        "--methods {input.methods:q} --systems {input.systems:q} "
        "--neutral-input {input.neutral_input:q} --neutral-output {input.neutral_output:q} "
        "--anion-input {input.anion_input:q} --anion-output {input.anion_output:q} "
        "--spin-cube {input.spin_cube:q} --electron-cube {input.electron_cube:q} "
        "--output {output:q}"


rule summarize_stage_b2c_method_benchmark:
    input:
        plan=rules.plan_stage_b2c_method_benchmark.output,
        records=STAGE_B2C_METHOD_BENCHMARK_ANALYSES,
        methods="configs/methods.yaml",
        script=STAGE_B2C_METHOD_BENCHMARK,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_method_benchmark.summary.json"
    params:
        production_summary=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential_production.summary.json"
        ),
        production_gate=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential_production.done"
        ),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} summary "
        "--campaign {CAMPAIGN:q} --production-summary {params.production_summary:q} "
        "--production-gate {params.production_gate:q} --plan {input.plan:q} "
        "--methods {input.methods:q} --records {input.records:q} --output {output:q}"


rule stage_b2c_method_benchmark:
    input:
        summary=rules.summarize_stage_b2c_method_benchmark.output,
        script=STAGE_B2C_METHOD_BENCHMARK,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_method_benchmark.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"
