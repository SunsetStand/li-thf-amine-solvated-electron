&GLOBAL
  PROJECT $project
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    CHARGE $charge
    MULTIPLICITY $multiplicity
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
    &KIND Gh
      ELEMENT H
      GHOST TRUE
      BASIS_SET $ghost_basis_set
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
