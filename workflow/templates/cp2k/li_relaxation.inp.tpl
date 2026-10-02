&GLOBAL
  PROJECT li_relax
  RUN_TYPE $run_type
  PRINT_LEVEL MEDIUM
&END GLOBAL
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    CHARGE 0
    MULTIPLICITY 2
    UKS TRUE
$basis_files
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF $cutoff_ry
      REL_CUTOFF $rel_cutoff_ry
      NGRIDS 5
    &END MGRID
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-12
      &CDFT
        TYPE_OF_CONSTRAINT BECKE
        ATOMIC_CHARGES TRUE
        STRENGTH 0.0
        TARGET $target
        &ATOM_GROUP
          ATOMS $li_atom_index
          COEFF 1.0
          CONSTRAINT_TYPE CHARGE
        &END ATOM_GROUP
        &OUTER_SCF ON
          TYPE CDFT_CONSTRAINT
          MAX_SCF 40
          EPS_SCF $cdft_eps_scf
          OPTIMIZER NEWTON_LS
          STEP_SIZE -1.0
          &CDFT_OPT ON
            MAX_LS 10
            CONTINUE_LS
            FACTOR_LS 0.5
            JACOBIAN_STEP 1.0E-2
            JACOBIAN_FREQ 1 1
            JACOBIAN_TYPE FD1
            JACOBIAN_RESTART FALSE
          &END CDFT_OPT
        &END OUTER_SCF
        &BECKE_CONSTRAINT
          CUTOFF_TYPE GLOBAL
          GLOBAL_CUTOFF 8.0
          CAVITY_CONFINE FALSE
          SHOULD_SKIP FALSE
        &END BECKE_CONSTRAINT
      &END CDFT
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF $eps_scf
      MAX_SCF $max_scf
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF ON
        EPS_SCF $eps_scf
        MAX_SCF 20
      &END OUTER_SCF
    &END SCF
$admm_block
    &XC
      &XC_FUNCTIONAL
        &PBE
          SCALE_X $pbe_scale_x
          SCALE_C 1.0
        &END PBE
      &END XC_FUNCTIONAL
$hf_block
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL $xc_reference
          R_CUTOFF 15.0
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
    &PRINT
$density_block
      &MULLIKEN ON
      &END MULLIKEN
    &END PRINT
  &END DFT
  &SUBSYS
    @INCLUDE $cell_path
    &TOPOLOGY
      COORD_FILE_NAME $coordinates_path
      COORD_FILE_FORMAT XYZ
      CONNECTIVITY OFF
    &END TOPOLOGY
$kinds
  &END SUBSYS
&END FORCE_EVAL
$motion_block
