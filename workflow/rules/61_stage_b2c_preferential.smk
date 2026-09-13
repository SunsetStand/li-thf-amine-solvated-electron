rule prepare_stage_b2c_preferential_seeds:
    input:
        methods="configs/methods.yaml",
        systems="configs/systems.yaml",
        script=PREPARE_STAGE_B2C_PREFERENTIAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        manifest=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/"
            "{system}/r{replica}/manifest.json"
        ),
        coordinates=[
            (
                f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/"
                f"{{system}}/r{{replica}}/{seed_role}/coordinates.xyz"
            )
            for seed_role in STAGE_B2C_SEED_ROLES
        ],
        cells=[
            (
                f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/"
                f"{{system}}/r{{replica}}/{seed_role}/cell.inc"
            )
            for seed_role in STAGE_B2C_SEED_ROLES
        ],
        metadata=[
            (
                f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/"
                f"{{system}}/r{{replica}}/{seed_role}/metadata.json"
            )
            for seed_role in STAGE_B2C_SEED_ROLES
        ],
    params:
        candidate_gate=stage_b2c_candidate_gate,
        source_campaign=lambda wildcards: stage_b2c_source_campaign(wildcards.system),
        candidate_summary=lambda wildcards: stage_b2c_artifact(
            f"{RUN_ROOT}/{stage_b2c_source_campaign(wildcards.system)}/"
            "stage_b_candidates.summary.json"
        ),
        candidate_manifest=lambda wildcards: stage_b2c_artifact(
            f"{RUN_ROOT}/{stage_b2c_source_campaign(wildcards.system)}/stage_b/"
            f"{wildcards.system}/r{wildcards.replica}/candidates/manifest.json"
        ),
        spec=lambda wildcards: stage_b2c_artifact(
            f"{RUN_ROOT}/{stage_b2c_source_campaign(wildcards.system)}/specs/"
            f"{wildcards.system}/r{wildcards.replica}.json"
        ),
        output_dir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/"
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
        "--campaign {CAMPAIGN:q} --source-campaign {params.source_campaign:q} "
        "--candidate-manifest {params.candidate_manifest:q} "
        "--candidate-summary {params.candidate_summary:q} "
        "--candidate-gate {params.candidate_gate:q} --spec {params.spec:q} "
        "--methods {input.methods:q} --systems {input.systems:q} "
        "--output-dir {params.output_dir:q} --output {output.manifest:q}"


rule render_stage_b2c_preferential_pair:
    input:
        manifest=rules.prepare_stage_b2c_preferential_seeds.output.manifest,
        coordinates=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/{wildcards.system}/"
            f"r{wildcards.replica}/{wildcards.seed_role}/coordinates.xyz"
        ),
        cell=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/{wildcards.system}/"
            f"r{wildcards.replica}/{wildcards.seed_role}/cell.inc"
        ),
        metadata=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/{wildcards.system}/"
            f"r{wildcards.replica}/{wildcards.seed_role}/metadata.json"
        ),
        methods="configs/methods.yaml",
        template="workflow/templates/cp2k/stage_b2c_preferential_smoke.inp.tpl",
        script=PREPARE_STAGE_B2C_PREFERENTIAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        neutral=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/"
            "{system}/r{replica}/{seed_role}/neutral/cp2k.inp"
        ),
        anion=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/"
            "{system}/r{replica}/{seed_role}/anion/cp2k.inp"
        ),
    params:
        project_prefix=lambda wildcards: (
            f"{wildcards.system}_r{wildcards.replica}_{wildcards.seed_role}"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed_role="eda_(rich|poor)",
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


rule run_stage_b2c_preferential_neutral:
    input:
        cp2k=rules.render_stage_b2c_preferential_pair.output.neutral,
        engine_runner=ENGINE_RUNNER,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/"
            "{system}/r{replica}/{seed_role}/neutral/cp2k.out"
        )
    params:
        workdir=lambda wildcards: stage_b2c_directory(
            wildcards.system, wildcards.replica, wildcards.seed_role, "neutral"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed_role="eda_(rich|poor)",
    threads: 1
    resources:
        tasks=32,
        mpi="mpirun",
        mem_mb=128000,
        runtime=720,
    shell:
        "bash {STAGE_RUNNER:q} cdft -- {PYTHON} {input.engine_runner:q} "
        "--engine cp2k --output {output:q} --cwd {params.workdir:q} -- "
        "{resources.mpi} -n {resources.tasks} cp2k.psmp -i {input.cp2k:q}"


rule run_stage_b2c_preferential_anion:
    input:
        cp2k=rules.render_stage_b2c_preferential_pair.output.anion,
        engine_runner=ENGINE_RUNNER,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        cp2k=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/"
            "{system}/r{replica}/{seed_role}/anion/cp2k.out"
        ),
        electron_cube=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/{{system}}/r{{replica}}/"
            "{seed_role}/anion/{system}_r{replica}_{seed_role}_anion_"
            "stage_b2c_preferential_smoke-ELECTRON_DENSITY-1_0.cube"
        ),
        spin_cube=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/{{system}}/r{{replica}}/"
            "{seed_role}/anion/{system}_r{replica}_{seed_role}_anion_"
            "stage_b2c_preferential_smoke-SPIN_DENSITY-1_0.cube"
        ),
    params:
        workdir=lambda wildcards: stage_b2c_directory(
            wildcards.system, wildcards.replica, wildcards.seed_role, "anion"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed_role="eda_(rich|poor)",
    threads: 1
    resources:
        tasks=32,
        mpi="mpirun",
        mem_mb=128000,
        runtime=720,
    shell:
        "bash {STAGE_RUNNER:q} cdft -- {PYTHON} {input.engine_runner:q} "
        "--engine cp2k --output {output.cp2k:q} --cwd {params.workdir:q} -- "
        "{resources.mpi} -n {resources.tasks} cp2k.psmp -i {input.cp2k:q}"


rule analyze_stage_b2c_preferential_pair:
    input:
        manifest=rules.prepare_stage_b2c_preferential_seeds.output.manifest,
        neutral_input=rules.render_stage_b2c_preferential_pair.output.neutral,
        anion_input=rules.render_stage_b2c_preferential_pair.output.anion,
        neutral_output=rules.run_stage_b2c_preferential_neutral.output,
        anion_output=rules.run_stage_b2c_preferential_anion.output.cp2k,
        electron_cube=rules.run_stage_b2c_preferential_anion.output.electron_cube,
        spin_cube=rules.run_stage_b2c_preferential_anion.output.spin_cube,
        methods="configs/methods.yaml",
        systems="configs/systems.yaml",
        script=PREPARE_STAGE_B2C_PREFERENTIAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential/"
            "{system}/r{replica}/{seed_role}/analysis.json"
        )
    params:
        spec=lambda wildcards: stage_b2c_artifact(
            f"{RUN_ROOT}/{stage_b2c_source_campaign(wildcards.system)}/specs/"
            f"{wildcards.system}/r{wildcards.replica}.json"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed_role="eda_(rich|poor)",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} analyze "
        "--campaign {CAMPAIGN:q} --system {wildcards.system:q} "
        "--replica {wildcards.replica} --seed-role {wildcards.seed_role:q} "
        "--manifest {input.manifest:q} --spec {params.spec:q} "
        "--methods {input.methods:q} --systems {input.systems:q} "
        "--neutral-input {input.neutral_input:q} --neutral-output {input.neutral_output:q} "
        "--anion-input {input.anion_input:q} --anion-output {input.anion_output:q} "
        "--spin-cube {input.spin_cube:q} --electron-cube {input.electron_cube:q} "
        "--output {output:q}"


rule summarize_stage_b2c_preferential_smoke:
    input:
        records=STAGE_B2C_ANALYSES,
        methods="configs/methods.yaml",
        script=PREPARE_STAGE_B2C_PREFERENTIAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential_smoke.summary.json"
    params:
        systems=" ".join(shlex.quote(value) for value in STAGE_B2C_SYSTEMS),
        replicas=" ".join(str(value) for value in STAGE_B2C_REPLICAS),
        seed_roles=" ".join(shlex.quote(value) for value in STAGE_B2C_SEED_ROLES),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} summary "
        "--campaign {CAMPAIGN:q} --methods {input.methods:q} --records {input.records:q} "
        "--expected-systems {params.systems} --expected-replicas {params.replicas} "
        "--expected-seed-roles {params.seed_roles} --output {output:q}"


rule stage_b2c_preferential_smoke:
    input:
        summary=rules.summarize_stage_b2c_preferential_smoke.output,
        script=PREPARE_STAGE_B2C_PREFERENTIAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b2c_preferential_smoke.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"
