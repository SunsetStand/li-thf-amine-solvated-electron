STAGE_B_REPORT_DIR = f"{RUN_ROOT}/{CAMPAIGN}/stage_b_report"


rule build_stage_b_report:
    input:
        script=BUILD_STAGE_B_REPORT,
        sources=SOLVELEC_SOURCES,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        pdf=f"{STAGE_B_REPORT_DIR}/stage_b_report_zh.pdf",
        metrics=f"{STAGE_B_REPORT_DIR}/stage_b_metrics.json",
        provenance=f"{STAGE_B_REPORT_DIR}/report_provenance.json",
        pipeline=f"{STAGE_B_REPORT_DIR}/figures/stage_b_pipeline.png",
        candidates=f"{STAGE_B_REPORT_DIR}/figures/candidate_bank.png",
        mechanism=f"{STAGE_B_REPORT_DIR}/figures/mechanism_energy.png",
        spin=f"{STAGE_B_REPORT_DIR}/figures/spin_partition.png",
        cubes=f"{STAGE_B_REPORT_DIR}/figures/cube_maps.png",
    params:
        campaign=CAMPAIGN,
        output_dir=STAGE_B_REPORT_DIR,
        # The completed Stage-B artifacts are intentionally an immutable handoff.
        # Keeping the gate in params prevents this report-only target from
        # traversing CP2K, GROMACS, or earlier analysis producers.
        accepted_stage_b=stage_b_report_handoff,
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} "
        "--run-root {RUN_ROOT:q} --campaign {params.campaign:q} "
        "--output-dir {params.output_dir:q}"


rule stage_b_report:
    input:
        report=rules.build_stage_b_report.output.pdf,
        metrics=rules.build_stage_b_report.output.metrics,
        provenance=rules.build_stage_b_report.output.provenance,
        script=PREPARE_STAGE_B,
        runtime=STAGE_RUNTIME_INPUTS,
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/stage_b_report.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bash {STAGE_RUNNER:q} trajectory_analysis -- {PYTHON} {input.script:q} gate "
        "--summary {input.provenance:q} --output {output:q}"
