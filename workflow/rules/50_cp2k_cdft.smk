rule render_cp2k_cdft:
    input:
        spec=f"{RUN_ROOT}/{CAMPAIGN}/specs/{{system}}/r{{replica}}.json",
        template="workflow/templates/cp2k/pbe0_cdft.inp.tpl"
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/inputs/{{system}}/r{{replica}}/cp2k_cdft.inp"
    params:
        project=lambda wildcards: f"{wildcards.system}_r{wildcards.replica}_cdft",
        xyz=lambda wildcards: f"qm/{wildcards.system}_r{wildcards.replica}.xyz",
        cell=lambda wildcards: f"qm/{wildcards.system}_r{wildcards.replica}.cell.inc"
    shell:
        "{PYTHON} -m solvelec.cli render-cp2k --state solvated_electron --constrained "
        "--project {params.project} --coordinates-include {params.xyz} "
        "--cell-include {params.cell} --li-atom-index 1 --output {output}"
