rule render_packmol:
    input:
        f"{RUN_ROOT}/{CAMPAIGN}/specs/{{system}}/r{{replica}}.json"
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/inputs/{{system}}/r{{replica}}/packmol.inp"
    params:
        packed=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/"
            f"r{wildcards.replica}/packed.pdb"
        ),
        thf=lambda wildcards: molecule_path("thf", "pdb"),
        amine_arg=packmol_amine_argument
    shell:
        "{PYTHON} -m solvelec.cli render-packmol --system {wildcards.system} "
        "--replica {wildcards.replica} --output {output:q} "
        "--output-pdb {params.packed:q} --thf-structure {params.thf:q} "
        "{params.amine_arg}"
