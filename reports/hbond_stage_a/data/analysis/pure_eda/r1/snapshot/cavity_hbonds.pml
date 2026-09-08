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
set sphere_scale, 1.9985, cavity_center
set sphere_transparency, 0.65, cavity_center
color cyan, cavity_center
distance hbond_001, shell and id 66, shell and id 88
color orange, hbond_001
set dash_width, 2, hbond_001
distance hbond_002, shell and id 83, shell and id 124
color orange, hbond_002
set dash_width, 2, hbond_002
distance hbond_003, shell and id 119, shell and id 85
color orange, hbond_003
set dash_width, 2, hbond_003
distance hbond_004, shell and id 137, shell and id 97
color orange, hbond_004
set dash_width, 2, hbond_004
distance hbond_005, shell and id 144, shell and id 148
color orange, hbond_005
set dash_width, 2, hbond_005
distance hbond_006, shell and id 180, shell and id 64
color orange, hbond_006
set dash_width, 2, hbond_006
set dash_gap, 0.18
set dash_length, 0.25
bg_color white
orient shell
zoom shell, 2.0
