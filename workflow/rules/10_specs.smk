rule composition_spec:
    input:
        rules.validate_config.output
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/specs/{{system}}/r{{replica}}.json"
    shell:
        "{PYTHON} -m solvelec.cli write-spec --system {wildcards.system} "
        "--replica {wildcards.replica} --output {output}"
