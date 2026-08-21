&GLOBAL
  PROJECT $project
  RUN_TYPE ENERGY_FORCE
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    CHARGE $charge
    MULTIPLICITY $multiplicity
    UKS TRUE
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    BASIS_SET_FILE_NAME BASIS_ADMM_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF $cutoff_ry
      REL_CUTOFF $rel_cutoff_ry
      NGRIDS 5
    &END MGRID
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-12
$cdft_block
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 100
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-6
        MAX_SCF 20
      &END OUTER_SCF
    &END SCF
    &AUXILIARY_DENSITY_MATRIX_METHOD
      ADMM_TYPE ADMMS
      EXCH_CORRECTION_FUNC PBEX
    &END AUXILIARY_DENSITY_MATRIX_METHOD
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &HF
        FRACTION $exact_exchange_fraction
        &SCREENING
          EPS_SCHWARZ 1.0E-6
          SCREEN_ON_INITIAL_P TRUE
        &END SCREENING
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE TRUNCATED
          CUTOFF_RADIUS $hfx_cutoff_angstrom
          T_C_G_DATA t_c_g.dat
        &END INTERACTION_POTENTIAL
        &MEMORY
          MAX_MEMORY 3000
          EPS_STORAGE_SCALING 0.1
        &END MEMORY
      &END HF
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE0
          R_CUTOFF 15.0
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
    &PRINT
      &E_DENSITY_CUBE
        STRIDE 1 1 1
      &END E_DENSITY_CUBE
      &SPIN_DENSITY_CUBE
        STRIDE 1 1 1
      &END SPIN_DENSITY_CUBE
      &HIRSHFELD ON
      &END HIRSHFELD
    &END PRINT
  &END DFT
  &SUBSYS
    @INCLUDE $cell_include
    &TOPOLOGY
      COORD_FILE_NAME $coordinates_include
      COORD_FILE_FORMAT XYZ
      CONNECTIVITY OFF
    &END TOPOLOGY
    &KIND H
      ELEMENT H
      BASIS_SET $basis_set
      BASIS_SET AUX_FIT $aux_basis_set
      POTENTIAL $potential
    &END KIND
    &KIND C
      ELEMENT C
      BASIS_SET $basis_set
      BASIS_SET AUX_FIT $aux_basis_set
      POTENTIAL $potential
    &END KIND
    &KIND N
      ELEMENT N
      BASIS_SET $basis_set
      BASIS_SET AUX_FIT $aux_basis_set
      POTENTIAL $potential
    &END KIND
    &KIND O
      ELEMENT O
      BASIS_SET $basis_set
      BASIS_SET AUX_FIT $aux_basis_set
      POTENTIAL $potential
    &END KIND
    &KIND Li
      ELEMENT Li
      BASIS_SET $basis_set
      BASIS_SET AUX_FIT $aux_basis_set
      POTENTIAL $potential
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
