#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
write_mdp()
write_continue_script()
make_gel(name='vacuum',**state.melt_settings)
minimize('vacuum')
solvate(structure='vacuum-minimized',gro='solvate')
write_topology('solvate.top')
copy_file('solvate.gro','system-input.gro')
copy_file('vacuum.top','system.top')
with open(state.here+'system.top','a') as fp:
	fp.write('%s %d\n'%(state.get('sol','SOL'),
		component(state.get('sol','SOL'))))
copy_file('system.top','solvate.top')
minimize('solvate')
copy_file('solvate-minimized.gro','system.gro')
equilibrate()
