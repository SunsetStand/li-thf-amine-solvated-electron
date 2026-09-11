rule analyze_stage_b_localization_state:
    input:
        methods="configs/methods.yaml",
        systems="configs/systems.yaml",
        script=ANALYZE_STAGE_B_LOCALIZATION,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        (
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{{system}}/r{{replica}}/mechanism_smoke/"
            "{candidate}/{state}/localization.json"
        )
    params:
        campaign=CAMPAIGN,
        # These files belong to the already accepted Stage-B mechanism run.
        # Keeping them in params prevents this analysis-only target from
        # traversing their producer rules and rebuilding prior calculations.
        mechanism_gate=stage_b_localization_mechanism_gate,
        spin_cube=lambda wildcards: stage_b_localization_artifact(
            stage_b_mechanism_cube(
                wildcards.system, wildcards.replica, wildcards.state, "SPIN_DENSITY"
            )
        ),
        electron_cube=lambda wildcards: stage_b_localization_artifact(
            stage_b_mechanism_cube(
                wildcards.system, wildcards.replica, wildcards.state, "ELECTRON_DENSITY"
            )
        ),
        candidate_metadata=lambda wildcards: stage_b_localization_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{wildcards.system}/r{wildcards.replica}/"
            f"candidates/{wildcards.candidate}/metadata.json"
        ),
        coordinates=lambda wildcards: stage_b_localization_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{wildcards.system}/r{wildcards.replica}/"
            f"candidates/{wildcards.candidate}/coordinates.xyz"
        ),
        cell=lambda wildcards: stage_b_localization_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{wildcards.system}/r{wildcards.replica}/"
            f"candidates/{wildcards.candidate}/cell.inc"
        ),
        spec=lambda wildcards: stage_b_localization_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/specs/{wildcards.system}/r{wildcards.replica}.json"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        candidate="[a-z][a-z0-9_-]*",
        state="[a-z][a-z0-9_]*",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} state "
        "--campaign {params.campaign:q} --system {wildcards.system:q} "
        "--replica {wildcards.replica} --candidate {wildcards.candidate:q} "
        "--state {wildcards.state:q} --spin-cube {params.spin_cube:q} "
        "--electron-cube {params.electron_cube:q} "
        "--candidate-metadata {params.candidate_metadata:q} "
        "--coordinates {params.coordinates:q} --cell {params.cell:q} --spec {params.spec:q} "
        "--methods {input.methods:q} --systems {input.systems:q} --output {output:q}"


rule analyze_stage_b_localization_pair:
    input:
        state_records=lambda wildcards: [
            (
                f"{stage_b_mechanism_directory(wildcards.system, wildcards.replica, state)}/"
                "localization.json"
            )
            for state in STAGE_B_MECHANISM_TARGETS
        ],
        methods="configs/methods.yaml",
        systems="configs/systems.yaml",
        script=ANALYZE_STAGE_B_LOCALIZATION,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        record=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{{system}}/r{{replica}}/mechanism_smoke/"
            "{candidate}/localization_pair.json"
        ),
        difference_cube=(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{{system}}/r{{replica}}/mechanism_smoke/"
            "{candidate}/electron_density_difference.cube"
        ),
    params:
        campaign=CAMPAIGN,
        mechanism_gate=stage_b_localization_mechanism_gate,
        li0_electron_cube=lambda wildcards: stage_b_localization_artifact(
            stage_b_mechanism_cube(
                wildcards.system, wildcards.replica, "li0_diabatic", "ELECTRON_DENSITY"
            )
        ),
        li_plus_e_electron_cube=lambda wildcards: stage_b_localization_artifact(
            stage_b_mechanism_cube(
                wildcards.system,
                wildcards.replica,
                "li_plus_e_diabatic",
                "ELECTRON_DENSITY",
            )
        ),
        candidate_metadata=lambda wildcards: stage_b_localization_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{wildcards.system}/r{wildcards.replica}/"
            f"candidates/{wildcards.candidate}/metadata.json"
        ),
        coordinates=lambda wildcards: stage_b_localization_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{wildcards.system}/r{wildcards.replica}/"
            f"candidates/{wildcards.candidate}/coordinates.xyz"
        ),
        cell=lambda wildcards: stage_b_localization_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b/{wildcards.system}/r{wildcards.replica}/"
            f"candidates/{wildcards.candidate}/cell.inc"
        ),
        spec=lambda wildcards: stage_b_localization_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/specs/{wildcards.system}/r{wildcards.replica}.json"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
        candidate="[a-z][a-z0-9_-]*",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} pair "
        "--campaign {params.campaign:q} --system {wildcards.system:q} "
        "--replica {wildcards.replica} --candidate {wildcards.candidate:q} "
        "--state-records {input.state_records:q} "
        "--li0-electron-cube {params.li0_electron_cube:q} "
        "--li-plus-e-electron-cube {params.li_plus_e_electron_cube:q} "
        "--difference-cube {output.difference_cube:q} "
        "--candidate-metadata {params.candidate_metadata:q} "
        "--coordinates {params.coordinates:q} --cell {params.cell:q} --spec {params.spec:q} "
        "--methods {input.methods:q} --systems {input.systems:q} --output {output.record:q}"


rule summarize_stage_b_localization:
    input:
        states=STAGE_B_LOCALIZATION_STATE_RECORDS,
        pairs=STAGE_B_LOCALIZATION_PAIR_RECORDS,
        script=ANALYZE_STAGE_B_LOCALIZATION,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        summary=f"{RUN_ROOT}/{CAMPAIGN}/stage_b_localization.summary.json",
        states_csv=f"{RUN_ROOT}/{CAMPAIGN}/stage_b_localization.states.csv",
        molecules_csv=f"{RUN_ROOT}/{CAMPAIGN}/stage_b_localization.molecules.csv",
        pairs_csv=f"{RUN_ROOT}/{CAMPAIGN}/stage_b_localization.pairs.csv",
    params:
        campaign=CAMPAIGN,
        # Delay this immutable-handoff check until the localization target is
        # actually selected.  Eager evaluation breaks unrelated smoke DAGs.
        mechanism_summary=lambda _wildcards: stage_b_localization_artifact(
            f"{RUN_ROOT}/{CAMPAIGN}/stage_b_mechanism_smoke.summary.json"
        ),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} summary "
        "--campaign {params.campaign:q} --mechanism-summary {params.mechanism_summary:q} "
        "--states {input.states:q} --pairs {input.pairs:q} "
        "--states-csv {output.states_csv:q} --molecules-csv {output.molecules_csv:q} "
        "--pairs-csv {output.pairs_csv:q} --output {output.summary:q}"


rule stage_b_localization:
    input:
        summary=rules.summarize_stage_b_localization.output.summary,
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b_localization.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"
