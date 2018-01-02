#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
forward_mapper()

#---begin parts from script-cg-gel on v027 which used the dextran_scale_switching/dextran_coarse_basic expt
write_mdp()
write_continue_script()
make_cg_gel(name='vacuum',**state.melt_settings)
write_topology('vacuum.top')
minimize('vacuum')
solvate(structure='vacuum-minimized',gro='solvate')
write_topology('solvate.top')
minimize('solvate')

copy_file('solvate.gro','system-input.gro')
#---! hacked
copy_file('vacuum.top','system.top')
with open(state.here+'system.top','a') as fp:
	fp.write('%s %d\n'%(state.q('sol','SOL'),component(state.q('sol','SOL'))))
copy_file('system.top','solvate.top')
#---! end hack
copy_file('solvate-minimized.gro','system.gro')
equilibrate()
