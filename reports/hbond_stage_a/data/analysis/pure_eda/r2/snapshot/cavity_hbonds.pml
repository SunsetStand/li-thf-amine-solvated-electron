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
set sphere_scale, 1.7955, cavity_center
set sphere_transparency, 0.65, cavity_center
color cyan, cavity_center
distance hbond_001, shell and id 6, shell and id 196
color orange, hbond_001
set dash_width, 2, hbond_001
distance hbond_002, shell and id 53, shell and id 181
color orange, hbond_002
set dash_width, 2, hbond_002
distance hbond_003, shell and id 77, shell and id 25
color orange, hbond_003
set dash_width, 2, hbond_003
distance hbond_004, shell and id 89, shell and id 145
color orange, hbond_004
set dash_width, 2, hbond_004
distance hbond_005, shell and id 102, shell and id 85
color orange, hbond_005
set dash_width, 2, hbond_005
distance hbond_006, shell and id 120, shell and id 184
color orange, hbond_006
set dash_width, 2, hbond_006
distance hbond_007, shell and id 125, shell and id 148
color orange, hbond_007
set dash_width, 2, hbond_007
distance hbond_008, shell and id 126, shell and id 76
color orange, hbond_008
set dash_width, 2, hbond_008
distance hbond_009, shell and id 132, shell and id 76
color orange, hbond_009
set dash_width, 2, hbond_009
distance hbond_010, shell and id 143, shell and id 64
color orange, hbond_010
set dash_width, 2, hbond_010
set dash_gap, 0.18
set dash_length, 0.25
bg_color white
orient shell
zoom shell, 2.0
