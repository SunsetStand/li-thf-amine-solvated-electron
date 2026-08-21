rule campaign_report:
    input:
        SPEC_OUTPUTS,
        rules.provenance_manifest.output
    output:
        f"runs/{CAMPAIGN}/report/README.md"
    shell:
        "{PYTHON} -m solvelec.cli report --campaign {CAMPAIGN} --output {output}"


rule input_bundle:
    input:
        PACKMOL_OUTPUTS,
        CP2K_OUTPUTS,
        rules.campaign_report.output
    output:
        touch(f"runs/{CAMPAIGN}/input_bundle.done")
