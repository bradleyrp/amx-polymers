display projection Orthographic
mol new s02-coarse/system_chains.gro
animate delete beg 0 end 0 skip 0 0
mol addfile s02-coarse/md-short1.xtc
mol addfile s02-coarse/md-short2.xtc
mol addfile s02-coarse/md-short3.xtc
mol addfile s02-coarse/md.part0001.xtc
source ~/libs/cg_bonds.tcl
cg_bonds -tpr s02-coarse/md-short1.tpr -gmx gmx -cutoff 10
mol delrep 0 top
mol selection "not resname W"
mol addrep top
mol modcolor 0 0 ResID

mol addrep top
mol modstyle 1 0 Points 6.000000
mol modcolor 1 0 ResID

# mol smoothrep 0 0 5
# mol smoothrep 0 1 5