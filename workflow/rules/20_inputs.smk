rule render_packmol:
    input:
        f"runs/{CAMPAIGN}/specs/{{system}}/r{{replica}}.json"
    output:
        f"runs/{CAMPAIGN}/inputs/{{system}}/r{{replica}}/packmol.inp"
    shell:
        "{PYTHON} -m solvelec.cli render-packmol --system {wildcards.system} "
        "--replica {wildcards.replica} --output {output}"
