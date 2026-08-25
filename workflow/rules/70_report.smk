rule campaign_report:
    input:
        SPEC_OUTPUTS,
        rules.provenance_manifest.output
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/report/README.md"
    shell:
        "{PYTHON} -m solvelec.cli report --campaign {CAMPAIGN} --output {output}"


rule input_bundle:
    input:
        PACKMOL_OUTPUTS,
        CP2K_OUTPUTS,
        rules.campaign_report.output
    output:
        touch(f"{RUN_ROOT}/{CAMPAIGN}/input_bundle.done")
