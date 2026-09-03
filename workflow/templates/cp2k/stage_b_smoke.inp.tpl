&GLOBAL
  PROJECT $project
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    CHARGE 0
    MULTIPLICITY 2
    UKS TRUE
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF $cutoff_ry
      REL_CUTOFF $rel_cutoff_ry
      NGRIDS 5
    &END MGRID
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-10
      &CDFT
        TYPE_OF_CONSTRAINT BECKE
        ATOMIC_CHARGES TRUE
        STRENGTH 0.0
        # TARGET is a valence-electron population. The installed Li
        # GTH-PBE-q3 potential therefore represents Li+ at 2 electrons.
        TARGET $li_target_valence_electrons
        &ATOM_GROUP
          ATOMS $li_atom_index
          COEFF 1.0
          CONSTRAINT_TYPE CHARGE
        &END ATOM_GROUP
        &OUTER_SCF ON
          TYPE CDFT_CONSTRAINT
          EXTRAPOLATION_ORDER 2
          MAX_SCF 20
          EPS_SCF $cdft_eps_scf
          OPTIMIZER NEWTON_LS
          STEP_SIZE -1.0
          &CDFT_OPT ON
            MAX_LS 5
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
          GLOBAL_CUTOFF 6.0
          CAVITY_CONFINE TRUE
          CAVITY_SHAPE VDW
          EPS_CAVITY 1.0E-7
          SHOULD_SKIP TRUE
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
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
          R_CUTOFF 15.0
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
    &PRINT
      &E_DENSITY_CUBE
        STRIDE $cube_stride $cube_stride $cube_stride
      &END E_DENSITY_CUBE
      &HIRSHFELD ON
      &END HIRSHFELD
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
    &KIND H
      ELEMENT H
      BASIS_SET $basis_set
      POTENTIAL $potential
    &END KIND
    &KIND C
      ELEMENT C
      BASIS_SET $basis_set
      POTENTIAL $potential
    &END KIND
    &KIND N
      ELEMENT N
      BASIS_SET $basis_set
      POTENTIAL $potential
    &END KIND
    &KIND O
      ELEMENT O
      BASIS_SET $basis_set
      POTENTIAL $potential
    &END KIND
    &KIND Li
      ELEMENT Li
      BASIS_SET $basis_set
      POTENTIAL $li_potential
    &END KIND
    &KIND Gh
      ELEMENT H
      GHOST TRUE
      BASIS_SET $ghost_basis_set
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
