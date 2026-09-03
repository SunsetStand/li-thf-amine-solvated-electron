rule prepare_stage_b_candidates:
    input:
        methods="configs/methods.yaml",
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    params:
        # Stage A is an immutable handoff. These paths deliberately remain
        # params so Stage B can never reschedule or rewrite solvent trajectories
        # or the accepted snapshot bank.
        spec=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/specs/{wildcards.system}/r{wildcards.replica}.json"
        ),
        snapshot_metadata=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{wildcards.system}/r{wildcards.replica}/"
            "snapshot/metadata.json"
        ),
        xyz=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{wildcards.system}/r{wildcards.replica}/"
            "snapshot/representative.xyz"
        ),
        cell=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{wildcards.system}/r{wildcards.replica}/"
            "snapshot/representative.cell.inc"
        ),
        output_dir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{wildcards.system}/r{wildcards.replica}/candidates"
        ),
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{{system}}/r{{replica}}/candidates/manifest.json"
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} prepare "
        "--spec {params.spec:q} --snapshot-metadata {params.snapshot_metadata:q} "
        "--xyz {params.xyz:q} --cell {params.cell:q} --methods {input.methods:q} "
        "--output-dir {params.output_dir:q} --output {output:q}"


rule summarize_stage_b_candidates:
    input:
        manifests=STAGE_B_MANIFESTS,
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b_candidates.summary.json"
    params:
        campaign=CAMPAIGN,
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} summary "
        "--campaign {params.campaign:q} --output {output:q} {input.manifests:q}"


rule stage_b_candidates:
    input:
        summary=rules.summarize_stage_b_candidates.output,
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b_candidates.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"


rule render_stage_b_cp2k_smoke:
    input:
        candidate_gate=rules.stage_b_candidates.output,
        manifest=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{{system}}/r{{replica}}/candidates/manifest.json"
        ),
        methods="configs/methods.yaml",
        template="workflow/templates/cp2k/stage_b_smoke.inp.tpl",
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{{system}}/r{{replica}}/smoke/"
            "{candidate}/cp2k.inp"
        )
    params:
        project=lambda wildcards: (
            f"{wildcards.system}_r{wildcards.replica}_{wildcards.candidate}_bsmoke"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        candidate="[a-z][a-z0-9_-]*",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} render "
        "--manifest {input.manifest:q} --methods {input.methods:q} "
        "--template {input.template:q} --candidate {wildcards.candidate:q} "
        "--project {params.project:q} --output {output:q}"


rule run_stage_b_cp2k_smoke:
    input:
        cp2k=rules.render_stage_b_cp2k_smoke.output,
        engine_runner=ENGINE_RUNNER,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{{system}}/r{{replica}}/smoke/"
            "{candidate}/cp2k.out"
        )
    params:
        workdir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{wildcards.system}/r{wildcards.replica}/"
            f"smoke/{wildcards.candidate}"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        candidate="[a-z][a-z0-9_-]*",
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


rule summarize_stage_b_cp2k_smoke:
    input:
        outputs=STAGE_B_SMOKE_OUTPUTS,
        cp2k_inputs=STAGE_B_SMOKE_INPUTS,
        manifests=STAGE_B_SMOKE_MANIFESTS,
        methods="configs/methods.yaml",
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b_cp2k_smoke.summary.json"
    params:
        campaign=CAMPAIGN,
        candidate=STAGE_B_SMOKE_CANDIDATE,
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} "
        "smoke-summary --campaign {params.campaign:q} --candidate {params.candidate:q} "
        "--methods {input.methods:q} --output {output:q} --outputs {input.outputs:q} "
        "--cp2k-inputs {input.cp2k_inputs:q} --manifests {input.manifests:q}"


rule stage_b:
    input:
        summary=rules.summarize_stage_b_cp2k_smoke.output,
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"
