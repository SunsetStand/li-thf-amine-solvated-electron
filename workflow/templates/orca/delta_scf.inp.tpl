! UKS $functional $basis_set $grid $scf NoAutoStart

%scf
  MaxIter 500
end

%output
  Print[P_MOs] 1
  Print[P_AtCharges_M] 1
end

* xyz $charge $multiplicity
$coordinates_xyz
*
