rule analyze_classical_replica:
    input:
        methods="configs/methods.yaml",
        script=ANALYZE_CLASSICAL,
        runtime=STAGE_RUNTIME_INPUTS,
        campaign_handoff=hbond_analysis_handoff,
    params:
        campaign=require_pilot_campaign,
        # Stage A consumes a completed, validated pilot as an immutable data
        # product. Keeping these paths in params deliberately prevents a newer
        # analysis config/source file from scheduling the 20 ns MD producers.
        spec=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/specs/{wildcards.system}/r{wildcards.replica}.json"
        ),
        classical_validation=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/"
            "pilot/validation.json"
        ),
        tpr=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/"
            "pilot/production/production.tpr"
        ),
        trajectory=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/"
            "pilot/production/production.xtc"
        ),
    output:
        analysis=f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/analysis.json",
        timeseries=f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/timeseries.csv",
        rdf=f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/rdf.csv",
        hydrogen_bonds=(
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/hydrogen_bonds.csv"
        ),
    wildcard_constraints:
        system="[a-z0-9_]+",
        replica="[1-9][0-9]*",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=180,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} analyze "
        "--spec {params.spec:q} --methods {input.methods:q} "
        "--classical-validation {params.classical_validation:q} --tpr {params.tpr:q} "
        "--trajectory {params.trajectory:q} --timeseries {output.timeseries:q} "
        "--rdf {output.rdf:q} --hydrogen-bonds {output.hydrogen_bonds:q} "
        "--output {output.analysis:q}"


rule summarize_classical_analysis:
    input:
        CLASSICAL_ANALYSIS_RECORDS,
        script=ANALYZE_CLASSICAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/classical_analysis.summary.json"
    params:
        campaign=require_pilot_campaign,
        records=lambda _wildcards: " ".join(
            shlex.quote(path) for path in CLASSICAL_ANALYSIS_RECORDS
        ),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} summary "
        "--campaign {params.campaign:q} --kind analysis --output {output:q} "
        "{params.records}"


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
        script=ANALYZE_CLASSICAL,
        runtime=STAGE_RUNTIME_INPUTS,
    params:
        tpr=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/"
            "pilot/production/production.tpr"
        ),
        trajectory=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/"
            "pilot/production/production.xtc"
        ),
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
        cavity_hbonds=(
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/snapshot/"
            "cavity_hbonds.json"
        ),
        local_pdb=(
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/snapshot/"
            "cavity_local.pdb"
        ),
        pymol=(
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/snapshot/"
            "cavity_hbonds.pml"
        ),
        svg=(
            f"{RUN_ROOT}/{CAMPAIGN}/analysis/{{system}}/r{{replica}}/snapshot/"
            "cavity_hbonds.svg"
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
        "--methods {input.methods:q} --tpr {params.tpr:q} --trajectory {params.trajectory:q} "
        "--xyz {output.xyz:q} --cell {output.cell:q} "
        "--cavity-hbonds {output.cavity_hbonds:q} --local-pdb {output.local_pdb:q} "
        "--pymol {output.pymol:q} --svg {output.svg:q} --output {output.metadata:q}"


rule summarize_snapshot_bank:
    input:
        SNAPSHOT_RECORDS,
        script=ANALYZE_CLASSICAL,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/snapshot_bank.summary.json"
    params:
        campaign=require_pilot_campaign,
        records=lambda _wildcards: " ".join(shlex.quote(path) for path in SNAPSHOT_RECORDS),
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} summary "
        "--campaign {params.campaign:q} --kind snapshot --output {output:q} "
        "{params.records}"


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


rule hbond_stage_a:
    input:
        classical=rules.classical_pilot.output,
        analysis=rules.classical_analysis.output,
        snapshots=rules.snapshot_bank.output,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/hbond_stage_a.done"
    params:
        campaign=require_hbond_pilot_campaign,
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "touch {output:q}"
