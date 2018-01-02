display projection Orthographic
mol new /home/rpb/omicron/dataset-project-polymers/melt-v011/s01-melt/system.gro
animate delete beg 0 end 0 skip 0 0
mol addfile /home/rpb/omicron/dataset-project-polymers/melt-v011/s01-melt/md.parts.pbcmol.centered.xtc step 100 waitfor all
mol delrep 0 top
mol selection "resname AGLC"
mol addrep top
mol new s01-tune-fine-to-coarse/cg_model.gro
mol modstyle 0 1 Points 16.000000
