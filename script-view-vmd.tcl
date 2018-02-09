display projection Orthographic
mol new s02-coarse/system.gro
animate delete beg 0 end 0 skip 0 0
mol addfile s02-coarse/md-short1.xtc
mol addfile s02-coarse/md-short2.xtc
mol addfile s02-coarse/md-short3.xtc
source ~/libs/cg_bonds.tcl
cg_bonds -tpr s02-coarse/md-short1.tpr -gmx gmx -cutoff 20
mol delrep 0 top
mol selection "not resname W"
mol addrep top
