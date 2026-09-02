rule analyze_classical_replica:
    input:
        spec=f"{RUN_ROOT}/{CAMPAIGN}/specs/{{system}}/r{{replica}}.json",
        methods="configs/methods.yaml",
        classical_validation=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/validation.json"
        ),
        tpr=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "production/production.tpr"
        ),
        trajectory=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "production/production.xtc"
        ),
        script=ANALYZE_CLASSICAL,
        runtime=STAGE_RUNTIME_INPUTS,
    params:
        campaign=require_pilot_campaign,
    output:
        analysis=f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/analysis.json",
        timeseries=f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/timeseries.csv",
        rdf=f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/rdf.csv",
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=180,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} analyze "
        "--spec {input.spec:q} --methods {input.methods:q} "
        "--classical-validation {input.classical_validation:q} --tpr {input.tpr:q} "
        "--trajectory {input.trajectory:q} --timeseries {output.timeseries:q} "
        "--rdf {output.rdf:q} --output {output.analysis:q}"


rule summarize_classical_analysis:
    input:
        CLASSICAL_ANALYSIS_RECORDS,
        script=ANALYZE_CLASSICAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/classical_analysis.summary.json"
    params:
        campaign=require_pilot_campaign,
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} summary "
        "--campaign {params.campaign:q} --kind analysis --output {output:q} "
        "{input[0]:q} {input[1]:q} {input[2]:q} {input[3]:q} {input[4]:q} {input[5]:q}"


rule classical_analysis:
    input:
        summary=rules.summarize_classical_analysis.output,
        script=ANALYZE_CLASSICAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/classical_analysis.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"


rule select_classical_snapshot:
    input:
        analysis_gate=rules.classical_analysis.output,
        analysis=f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/analysis.json",
        timeseries=f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/timeseries.csv",
        methods="configs/methods.yaml",
        tpr=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "production/production.tpr"
        ),
        trajectory=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "production/production.xtc"
        ),
        script=ANALYZE_CLASSICAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        xyz=(
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/snapshot/"
            "representative.xyz"
        ),
        cell=(
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/snapshot/"
            "representative.cell.inc"
        ),
        metadata=(
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/snapshot/metadata.json"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} select "
        "--analysis {input.analysis:q} --timeseries {input.timeseries:q} "
        "--methods {input.methods:q} --tpr {input.tpr:q} --trajectory {input.trajectory:q} "
        "--xyz {output.xyz:q} --cell {output.cell:q} --output {output.metadata:q}"


rule summarize_snapshot_bank:
    input:
        SNAPSHOT_RECORDS,
        script=ANALYZE_CLASSICAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/snapshot_bank.summary.json"
    params:
        campaign=require_pilot_campaign,
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} summary "
        "--campaign {params.campaign:q} --kind snapshot --output {output:q} "
        "{input[0]:q} {input[1]:q} {input[2]:q} {input[3]:q} {input[4]:q} {input[5]:q}"


rule snapshot_bank:
    input:
        summary=rules.summarize_snapshot_bank.output,
        script=ANALYZE_CLASSICAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/snapshot_bank.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.summary:q} --output {output:q}"
