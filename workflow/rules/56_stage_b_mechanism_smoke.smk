rule render_stage_b_mechanism_smoke:
    input:
        candidate_gate=stage_b_mechanism_candidate_gate,
        manifest=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{{system}}/r{{replica}}/candidates/manifest.json"
        ),
        methods="configs/methods.yaml",
        template="workflow/templates/cp2k/stage_b_smoke.inp.tpl",
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{{system}}/r{{replica}}/mechanism_smoke/"
            "{candidate}/{state}/cp2k.inp"
        )
    params:
        project=lambda wildcards: (
            f"{wildcards.system}_r{wildcards.replica}_{wildcards.candidate}_"
            f"{wildcards.state}_mechanism_smoke"
        ),
        target=lambda wildcards: STAGE_B_MECHANISM_TARGETS[wildcards.state],
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        candidate="[a-z][a-z0-9_-]*",
        state="[a-z][a-z0-9_]*",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} render "
        "--manifest {input.manifest:q} --methods {input.methods:q} "
        "--template {input.template:q} --candidate {wildcards.candidate:q} "
        "--project {params.project:q} --target-electrons {params.target} --output {output:q}"


rule run_stage_b_mechanism_smoke:
    input:
        cp2k=rules.render_stage_b_mechanism_smoke.output,
        engine_runner=ENGINE_RUNNER,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{{system}}/r{{replica}}/mechanism_smoke/"
            "{candidate}/{state}/cp2k.out"
        )
    params:
        workdir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{wildcards.system}/r{wildcards.replica}/"
            f"mechanism_smoke/{wildcards.candidate}/{wildcards.state}"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        candidate="[a-z][a-z0-9_-]*",
        state="[a-z][a-z0-9_]*",
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


rule summarize_stage_b_mechanism_smoke:
    input:
        outputs=STAGE_B_MECHANISM_OUTPUTS,
        cp2k_inputs=STAGE_B_MECHANISM_INPUTS,
        manifests=STAGE_B_MECHANISM_MANIFESTS,
        methods="configs/methods.yaml",
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b_mechanism_smoke.summary.json"
    params:
        campaign=CAMPAIGN,
        candidate=STAGE_B_SMOKE_CANDIDATE,
        states=" ".join(shlex.quote(value) for value in STAGE_B_MECHANISM_STATE_IDS),
        targets=" ".join(str(value) for value in STAGE_B_MECHANISM_TARGET_VALUES),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} "
        "mechanism-summary --campaign {params.campaign:q} --candidate {params.candidate:q} "
        "--methods {input.methods:q} --output {output:q} --outputs {input.outputs:q} "
        "--cp2k-inputs {input.cp2k_inputs:q} --manifests {input.manifests:q} "
        "--states {params.states} --targets {params.targets}"


rule stage_b_mechanism_smoke:
    input:
        summary=rules.summarize_stage_b_mechanism_smoke.output,
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b_mechanism_smoke.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"
