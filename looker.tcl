mol new s01-melt/system.gro
mol addfile s01-melt/md-short.xtc
mol addfile s01-melt/md.part0001.xtc
display projection orthographic
mol delrep 0 0
mol color Name
mol representation Lines 1.000000
mol selection not resname W
mol material Opaque
mol addrep 0
source ~/libs/cg_bonds.tcl
cg_bonds -tpr s01-melt/md-short.tpr -gmx ~/libs/gmxdump
mol modstyle 0 0 Lines 2.000000
mol modcolor 0 0 ResID