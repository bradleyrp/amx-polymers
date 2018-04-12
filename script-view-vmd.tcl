display projection Orthographic
set spot s02-large
mol new $spot/system.gro
animate delete beg 0 end 0 skip 0 0
mol addfile $spot/md-short1.xtc
mol addfile $spot/md-short2.xtc
mol addfile $spot/md-short3.xtc
mol addfile $spot/md.part0001.xtc
source ~/libs/cg_bonds.tcl
cg_bonds -tpr $spot/md-short1.tpr -gmx gmx -cutoff 10
mol delrep 0 top
mol selection "not resname W"

mol addrep top
mol modstyle 0 0 Points 6.000000
mol modcolor 0 0 ResID

# mol addrep top
# mol modcolor 1 0 ResID

# mol smoothrep 0 0 5
# mol smoothrep 0 1 5
