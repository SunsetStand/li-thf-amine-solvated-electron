rule prepare_stage_b2_intrinsic_seeds:
    input:
        methods="configs/methods.yaml",
        script=PREPARE_STAGE_B2_INTRINSIC,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        manifest=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{{system}}/r{{replica}}/manifest.json"
        ),
        coordinates=[
            (
                f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{{system}}/r{{replica}}/"
                f"{seed}/coordinates.xyz"
            )
            for seed in STAGE_B2_INTRINSIC_SEEDS
        ],
        cells=[
            (
                f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{{system}}/r{{replica}}/"
                f"{seed}/cell.inc"
            )
            for seed in STAGE_B2_INTRINSIC_SEEDS
        ],
        metadata=[
            (
                f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{{system}}/r{{replica}}/"
                f"{seed}/metadata.json"
            )
            for seed in STAGE_B2_INTRINSIC_SEEDS
        ],
    params:
        candidate_gate=stage_b2_intrinsic_candidate_gate,
        candidate_summary=lambda _wildcards: stage_b2_intrinsic_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b_candidates.summary.json"
        ),
        candidate_manifest=lambda wildcards: stage_b2_intrinsic_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{wildcards.system}/r{wildcards.replica}/"
            "candidates/manifest.json"
        ),
        spec=lambda wildcards: stage_b2_intrinsic_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/specs/{wildcards.system}/r{wildcards.replica}.json"
        ),
        output_dir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{wildcards.system}/"
            f"r{wildcards.replica}"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} prepare "
        "--candidate-manifest {params.candidate_manifest:q} "
        "--candidate-summary {params.candidate_summary:q} "
        "--candidate-gate {params.candidate_gate:q} --spec {params.spec:q} "
        "--methods {input.methods:q} --output-dir {params.output_dir:q} "
        "--output {output.manifest:q}"


rule render_stage_b2_intrinsic_smoke:
    input:
        manifest=rules.prepare_stage_b2_intrinsic_seeds.output.manifest,
        coordinates=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{wildcards.system}/"
            f"r{wildcards.replica}/{wildcards.seed}/coordinates.xyz"
        ),
        cell=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{wildcards.system}/"
            f"r{wildcards.replica}/{wildcards.seed}/cell.inc"
        ),
        metadata=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{wildcards.system}/"
            f"r{wildcards.replica}/{wildcards.seed}/metadata.json"
        ),
        methods="configs/methods.yaml",
        template="workflow/templates/cp2k/stage_b2_intrinsic_smoke.inp.tpl",
        script=PREPARE_STAGE_B2_INTRINSIC,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{{system}}/r{{replica}}/"
            "{seed}/cp2k.inp"
        )
    params:
        project=lambda wildcards: stage_b2_intrinsic_project(
            wildcards.system, wildcards.replica, wildcards.seed
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed="void_[0-9][0-9]",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} render "
        "--manifest {input.manifest:q} --methods {input.methods:q} "
        "--template {input.template:q} --seed {wildcards.seed:q} "
        "--project {params.project:q} --output {output:q}"


rule run_stage_b2_intrinsic_smoke:
    input:
        cp2k=rules.render_stage_b2_intrinsic_smoke.output,
        engine_runner=ENGINE_RUNNER,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        cp2k=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{{system}}/r{{replica}}/"
            "{seed}/cp2k.out"
        ),
        electron_cube=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{{system}}/r{{replica}}/{{seed}}/"
            "{system}_r{replica}_{seed}_stage_b2_intrinsic_smoke-"
            "ELECTRON_DENSITY-1_0.cube"
        ),
        spin_cube=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{{system}}/r{{replica}}/{{seed}}/"
            "{system}_r{replica}_{seed}_stage_b2_intrinsic_smoke-SPIN_DENSITY-1_0.cube"
        ),
    params:
        workdir=lambda wildcards: stage_b2_intrinsic_directory(
            wildcards.system, wildcards.replica, wildcards.seed
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed="void_[0-9][0-9]",
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


rule analyze_stage_b2_intrinsic_state:
    input:
        cp2k_input=rules.render_stage_b2_intrinsic_smoke.output,
        cp2k_output=rules.run_stage_b2_intrinsic_smoke.output.cp2k,
        electron_cube=rules.run_stage_b2_intrinsic_smoke.output.electron_cube,
        spin_cube=rules.run_stage_b2_intrinsic_smoke.output.spin_cube,
        manifest=rules.prepare_stage_b2_intrinsic_seeds.output.manifest,
        methods="configs/methods.yaml",
        systems="configs/systems.yaml",
        script=PREPARE_STAGE_B2_INTRINSIC,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic/{{system}}/r{{replica}}/"
            "{seed}/analysis.json"
        )
    params:
        campaign=CAMPAIGN,
        spec=lambda wildcards: stage_b2_intrinsic_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/specs/{wildcards.system}/r{wildcards.replica}.json"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        seed="void_[0-9][0-9]",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} analyze "
        "--campaign {params.campaign:q} --system {wildcards.system:q} "
        "--replica {wildcards.replica} --seed {wildcards.seed:q} "
        "--manifest {input.manifest:q} --spec {params.spec:q} "
        "--methods {input.methods:q} --systems {input.systems:q} "
        "--cp2k-input {input.cp2k_input:q} --cp2k-output {input.cp2k_output:q} "
        "--spin-cube {input.spin_cube:q} --electron-cube {input.electron_cube:q} "
        "--output {output:q}"


rule summarize_stage_b2_intrinsic_smoke:
    input:
        records=STAGE_B2_INTRINSIC_ANALYSES,
        methods="configs/methods.yaml",
        script=PREPARE_STAGE_B2_INTRINSIC,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic_smoke.summary.json"
    params:
        campaign=CAMPAIGN,
        systems=" ".join(shlex.quote(value) for value in STAGE_B2_INTRINSIC_SYSTEMS),
        replicas=" ".join(str(value) for value in STAGE_B2_INTRINSIC_REPLICAS),
        seeds=" ".join(shlex.quote(value) for value in STAGE_B2_INTRINSIC_SEEDS),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} summary "
        "--campaign {params.campaign:q} --methods {input.methods:q} "
        "--records {input.records:q} --expected-systems {params.systems} "
        "--expected-replicas {params.replicas} --expected-seeds {params.seeds} "
        "--output {output:q}"


rule stage_b2_intrinsic_smoke:
    input:
        summary=rules.summarize_stage_b2_intrinsic_smoke.output,
        script=PREPARE_STAGE_B2_INTRINSIC,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b2_intrinsic_smoke.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"
