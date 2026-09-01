wildcard_constraints:
    molecule="thf|eda|12pda|13pda|deta|tmeda",
    system="[a-z0-9_]+",
    replica="[1-9][0-9]*",
    stage="em|nvt|npt|production"


rule prepare_molecule:
    input:
        catalog="configs/systems.yaml",
        script=PREPARE_MOLECULE,
        runtime=STAGE_RUNTIME_INPUTS
    output:
        mol2=f"{SHARED_ROOT}/molecules/{{molecule}}/{{molecule}}.mol2",
        frcmod=f"{SHARED_ROOT}/molecules/{{molecule}}/{{molecule}}.frcmod",
        pdb=f"{SHARED_ROOT}/molecules/{{molecule}}/{{molecule}}.pdb",
        sdf=f"{SHARED_ROOT}/molecules/{{molecule}}/{{molecule}}.sdf",
        manifest=f"{SHARED_ROOT}/molecules/{{molecule}}/manifest.json",
        run_log=f"{SHARED_ROOT}/molecules/{{molecule}}/parameterize.log"
    params:
        output_dir=lambda wildcards: f"{SHARED_ROOT}/molecules/{wildcards.molecule}",
        smiles=lambda wildcards: molecule_record(wildcards.molecule)["smiles"],
        residue=lambda wildcards: molecule_record(wildcards.molecule)["residue"],
        charge=lambda wildcards: molecule_record(wildcards.molecule)["charge"]
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120
    shell:
        "bash {STAGE_RUNNER:q} molecule_generation -- {AMBER_PYTHON:q} "
        "{PREPARE_MOLECULE:q} --name {wildcards.molecule:q} "
        "--smiles {params.smiles:q} --residue {params.residue:q} "
        "--charge {params.charge} --output-dir {params.output_dir:q}"


rule run_packmol:
    input:
        specification=f"{RUN_ROOT}/{CAMPAIGN}/specs/{{system}}/r{{replica}}.json",
        packmol=f"{RUN_ROOT}/{CAMPAIGN}/inputs/{{system}}/r{{replica}}/packmol.inp",
        molecules=molecule_pdb_inputs,
        runtime=STAGE_RUNTIME_INPUTS
    output:
        packed=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/packed.pdb",
        run_log=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/packmol.log",
        validation=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/packmol.validation.json"
        )
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        "bash {STAGE_RUNNER:q} molecule_generation -- packmol "
        "< {input.packmol:q} > {output.run_log:q} 2>&1 && "
        "test -s {output.packed:q} && "
        "{PYTHON} -m solvelec.cli parse-output --engine packmol {output.run_log:q} "
        "> {output.validation:q}"


rule build_gromacs_system:
    input:
        specification=f"{RUN_ROOT}/{CAMPAIGN}/specs/{{system}}/r{{replica}}.json",
        packed=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/packed.pdb",
        parameters=molecule_parameter_inputs,
        script=BUILD_GROMACS,
        runtime=STAGE_RUNTIME_INPUTS
    output:
        topology=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/build/topol.top",
        coordinates=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/build/conf.gro",
        prmtop=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/build/system.prmtop",
        inpcrd=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/build/system.inpcrd",
        manifest=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/build/manifest.json"
    log:
        f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/build/tleap.log"
    params:
        output_dir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/build"
        ),
        molecule_args=tleap_molecule_arguments
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120
    shell:
        "bash {STAGE_RUNNER:q} molecule_generation -- {AMBER_PYTHON:q} "
        "{BUILD_GROMACS:q} --packed-pdb {input.packed:q} "
        "--spec {input.specification:q} --output-dir {params.output_dir:q} "
        "{params.molecule_args}"


rule render_gromacs_smoke_mdp:
    input:
        template=lambda wildcards: f"workflow/templates/gromacs/smoke_{wildcards.stage}.mdp.tpl",
        sources=SOLVELEC_SOURCES
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/{{stage}}/{{stage}}.mdp"
    wildcard_constraints:
        stage="nvt|npt"
    shell:
        "{PYTHON} -m solvelec.cli render-gromacs-mdp --stage {wildcards.stage} "
        "--replica {wildcards.replica} --output {output:q}"


rule gromacs_smoke_em:
    input:
        coordinates=rules.build_gromacs_system.output.coordinates,
        topology=rules.build_gromacs_system.output.topology,
        mdp="workflow/templates/gromacs/smoke_em.mdp",
        script=RUN_GROMACS,
        runtime=STAGE_RUNTIME_INPUTS
    output:
        tpr=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/em/em.tpr",
        coordinates=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/em/em.gro",
        energy=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/em/em.edr",
        engine_log=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/em/em.log",
        grompp_log=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/em/grompp.log",
        stdout_log=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/em/mdrun.stdout.log"
        ),
        manifest=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/em/manifest.json",
        validation=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/em/validation.json"
        )
    params:
        output_dir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/em"
        )
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60
    shell:
        "bash {STAGE_RUNNER:q} classical_md -- {PYTHON} {RUN_GROMACS:q} --phase em "
        "--mdp {input.mdp:q} --coordinates {input.coordinates:q} "
        "--topology {input.topology:q} --output-dir {params.output_dir:q} && "
        "{PYTHON} -m solvelec.cli parse-output --engine gromacs {output.engine_log:q} "
        "> {output.validation:q}"


rule gromacs_smoke_nvt:
    input:
        coordinates=rules.gromacs_smoke_em.output.coordinates,
        topology=rules.build_gromacs_system.output.topology,
        mdp=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/nvt/nvt.mdp"
        ),
        script=RUN_GROMACS,
        runtime=STAGE_RUNTIME_INPUTS
    output:
        tpr=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/nvt/nvt.tpr",
        coordinates=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/nvt/nvt.gro",
        energy=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/nvt/nvt.edr",
        checkpoint=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/nvt/nvt.cpt",
        engine_log=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/nvt/nvt.log",
        grompp_log=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/nvt/grompp.log",
        stdout_log=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/nvt/mdrun.stdout.log"
        ),
        manifest=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/nvt/manifest.json",
        validation=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/nvt/validation.json"
        )
    params:
        output_dir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/nvt"
        )
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60
    shell:
        "bash {STAGE_RUNNER:q} classical_md -- {PYTHON} {RUN_GROMACS:q} --phase nvt "
        "--mdp {input.mdp:q} --coordinates {input.coordinates:q} "
        "--topology {input.topology:q} --output-dir {params.output_dir:q} && "
        "{PYTHON} -m solvelec.cli parse-output --engine gromacs {output.engine_log:q} "
        "> {output.validation:q}"


rule gromacs_smoke_npt:
    input:
        coordinates=rules.gromacs_smoke_nvt.output.coordinates,
        checkpoint=rules.gromacs_smoke_nvt.output.checkpoint,
        topology=rules.build_gromacs_system.output.topology,
        mdp=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/npt/npt.mdp"
        ),
        script=RUN_GROMACS,
        runtime=STAGE_RUNTIME_INPUTS
    output:
        tpr=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/npt/npt.tpr",
        coordinates=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/npt/npt.gro",
        energy=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/npt/npt.edr",
        checkpoint=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/npt/npt.cpt",
        engine_log=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/npt/npt.log",
        grompp_log=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/npt/grompp.log",
        stdout_log=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/npt/mdrun.stdout.log"
        ),
        manifest=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/npt/manifest.json",
        validation=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/npt/validation.json"
        )
    params:
        output_dir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/npt"
        )
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60
    shell:
        "bash {STAGE_RUNNER:q} classical_md -- {PYTHON} {RUN_GROMACS:q} --phase npt "
        "--mdp {input.mdp:q} --coordinates {input.coordinates:q} "
        "--topology {input.topology:q} --checkpoint {input.checkpoint:q} "
        "--output-dir {params.output_dir:q} && "
        "{PYTHON} -m solvelec.cli parse-output --engine gromacs {output.engine_log:q} "
        "> {output.validation:q}"


rule classical_smoke:
    input:
        CLASSICAL_SMOKE_OUTPUTS,
        CLASSICAL_SMOKE_VALIDATIONS
    output:
        touch(f"{RUN_ROOT}/{CAMPAIGN}/classical_smoke.done")


rule render_gromacs_pilot_mdp:
    input:
        template=lambda wildcards: (
            f"workflow/templates/gromacs/pilot_{wildcards.stage}.mdp.tpl"
        ),
        campaign="configs/campaign.yaml",
        methods="configs/methods.yaml",
        sources=SOLVELEC_SOURCES
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/{{stage}}/{{stage}}.mdp"
    wildcard_constraints:
        stage="nvt|npt|production"
    shell:
        "{PYTHON} -m solvelec.cli render-gromacs-mdp --protocol pilot "
        "--stage {wildcards.stage} --system {wildcards.system:q} "
        "--replica {wildcards.replica} --output {output:q}"


rule gromacs_pilot_em:
    input:
        coordinates=rules.build_gromacs_system.output.coordinates,
        topology=rules.build_gromacs_system.output.topology,
        mdp="workflow/templates/gromacs/smoke_em.mdp",
        script=RUN_GROMACS,
        runtime=STAGE_RUNTIME_INPUTS
    output:
        tpr=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/em/em.tpr",
        coordinates=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/em/em.gro"
        ),
        energy=f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/em/em.edr",
        engine_log=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/em/em.log"
        ),
        grompp_log=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/em/grompp.log"
        ),
        stdout_log=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/em/"
            "mdrun.stdout.log"
        ),
        manifest=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/em/manifest.json"
        )
    params:
        output_dir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/pilot/em"
        )
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120
    shell:
        "bash {STAGE_RUNNER:q} classical_md -- {PYTHON} {RUN_GROMACS:q} --phase em "
        "--mdp {input.mdp:q} --coordinates {input.coordinates:q} "
        "--topology {input.topology:q} --output-dir {params.output_dir:q}"


rule gromacs_pilot_md:
    input:
        coordinates=pilot_previous_coordinates,
        checkpoint=pilot_previous_checkpoint,
        previous_validation=pilot_previous_validation,
        topology=rules.build_gromacs_system.output.topology,
        mdp=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/{stage}.mdp"
        ),
        script=RUN_GROMACS,
        runtime=STAGE_RUNTIME_INPUTS
    output:
        tpr=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/{stage}.tpr"
        ),
        coordinates=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/{stage}.gro"
        ),
        energy=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/{stage}.edr"
        ),
        checkpoint=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/{stage}.cpt"
        ),
        trajectory=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/{stage}.xtc"
        ),
        engine_log=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/{stage}.log"
        ),
        grompp_log=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/grompp.log"
        ),
        stdout_log=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/mdrun.stdout.log"
        ),
        manifest=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/manifest.json"
        )
    params:
        output_dir=lambda wildcards: (
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{wildcards.system}/r{wildcards.replica}/"
            f"pilot/{wildcards.stage}"
        ),
        checkpoint_arg=pilot_checkpoint_argument,
        checkpoint_minutes=lambda _wildcards: int(
            json.loads(Path("configs/methods.yaml").read_text())["classical_md"][
                "checkpoint_minutes"
            ]
        )
    threads: 4
    resources:
        mem_mb=8000,
        runtime=720
    wildcard_constraints:
        stage="nvt|npt|production"
    shell:
        "bash {STAGE_RUNNER:q} classical_md -- {PYTHON} {RUN_GROMACS:q} "
        "--phase {wildcards.stage:q} --mdp {input.mdp:q} "
        "--coordinates {input.coordinates:q} --topology {input.topology:q} "
        "{params.checkpoint_arg} --checkpoint-minutes {params.checkpoint_minutes} "
        "--restartable --expect-trajectory --output-dir {params.output_dir:q}"


rule validate_gromacs_pilot_stage:
    input:
        engine_log=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/{stage}.log"
        ),
        coordinates=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "{stage}/{stage}.gro"
        )
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
        "{stage}/engine.validation.json"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        "{PYTHON} -m solvelec.cli parse-output --engine gromacs {input.engine_log:q} "
        "> {output:q}"


rule validate_classical_pilot_replica:
    input:
        spec=f"{RUN_ROOT}/{CAMPAIGN}/specs/{{system}}/r{{replica}}.json",
        methods="configs/methods.yaml",
        tpr=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "production/production.tpr"
        ),
        trajectory=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "production/production.xtc"
        ),
        engine_validation=(
            f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/"
            "production/engine.validation.json"
        ),
        script=VALIDATE_CLASSICAL
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/classical/{{system}}/r{{replica}}/pilot/validation.json"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60
    shell:
        "bash {STAGE_RUNNER:q} classical_md -- {PYTHON} {VALIDATE_CLASSICAL:q} replica "
        "--spec {input.spec:q} --methods {input.methods:q} --tpr {input.tpr:q} "
        "--trajectory {input.trajectory:q} --engine-validation {input.engine_validation:q} "
        "--concentration-tolerance-m {config[concentration_tolerance_m]} "
        "--output {output:q}"


rule summarize_classical_pilot:
    input:
        CLASSICAL_PILOT_VALIDATIONS
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/classical_pilot.validation.json"
    params:
        campaign=require_pilot_campaign,
        inputs=lambda _wildcards: " ".join(
            shlex.quote(path) for path in CLASSICAL_PILOT_VALIDATIONS
        )
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        "{PYTHON} {VALIDATE_CLASSICAL:q} summary --campaign {params.campaign:q} "
        "--methods configs/methods.yaml --output {output:q} {params.inputs}"


rule classical_pilot:
    input:
        summary=rules.summarize_classical_pilot.output
    output:
        f"{RUN_ROOT}/{CAMPAIGN}/classical_pilot.done"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60
    shell:
        "{PYTHON} {VALIDATE_CLASSICAL:q} gate --summary {input.summary:q} --output {output:q}"
