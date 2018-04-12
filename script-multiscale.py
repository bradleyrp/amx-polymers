#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
forward_mapper()
# make a new CG simulation (mimics the method from script-cg-gel.py for experiment dextran_coarse_basic)
#! eventually put this in a loop
state.forward_dn = str(state.here)
make_step('coarse')
#! hardcoded name
shutil.copyfile(state.forward_dn+'dextran.itp',state.here+'dextran.itp')
shutil.copyfile(state.forward_dn+'dextran_constraints.itp',state.here+'dextran_constraints.itp')
write_mdp()
write_continue_script()
if not state.q('lattice_melt_settings',False): 
	make_polymer(name='vacuum')
else: 
	make_crude_coarse_polymer(name='vacuum')
	state.itp = ['dextran.itp']
write_topology('vacuum.top')
minimize('vacuum')
solvate(structure='vacuum-minimized',gro='solvate-untrimmed')
hydration_adjust(structure='solvate-untrimmed',gro='solvate')
write_topology('solvate.top')
minimize('solvate')
copy_file('solvate.gro','system-input.gro')
copy_file('vacuum.top','system.top')
with open(state.here+'system.top','a') as fp:
	fp.write('%s %d\n'%(state.q('sol','SOL'),component(state.q('sol','SOL'))))
copy_file('system.top','solvate.top')
copy_file('solvate-minimized.gro','system-residues.gro')
if settings.get('do_constraints',False):
	state.itp = ['dextran_constraints.itp']
	write_topology('system.top')
write_structure_by_chain(structure='system-residues',gro='system')
equilibrate()

# rm -rf s02-large state.json && cp state_1.json state.json && cp expt_2.json expt.json && python script_2.py
