rule validate_config:
    input:
        "configs/campaign.yaml",
        "configs/systems.yaml",
        "configs/methods.yaml",
        SOLVELEC_SOURCES
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/meta/config.valid.json"
    shell:
        "{PYTHON} -m solvelec.cli validate > {output}"


rule provenance_manifest:
    input:
        rules.validate_config.output,
        "configs/campaign.yaml",
        "configs/systems.yaml",
        "configs/methods.yaml"
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/meta/manifest.json"
    shell:
        "{PYTHON} -m solvelec.cli manifest --campaign {CAMPAIGN} --output {output} {input}"
