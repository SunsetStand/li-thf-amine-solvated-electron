load cavity_local.pdb, shell
hide everything, all
show sticks, shell
set stick_radius, 0.12, shell
color gray70, elem C
color red, elem O
color marine, elem N
color white, elem H
pseudoatom cavity_center, pos=[0.0, 0.0, 0.0]
show spheres, cavity_center
set sphere_scale, 1.8162, cavity_center
set sphere_transparency, 0.65, cavity_center
color cyan, cavity_center
set dash_gap, 0.18
set dash_length, 0.25
bg_color white
orient shell
zoom shell, 2.0
