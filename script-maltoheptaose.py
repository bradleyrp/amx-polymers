#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
write_mdp()
make_single_polymer()
component(settings.molecule_name,count=1)
state.itp = ['martini_v2.0_sugars.itp']
write_topology('vacuum_crude.top')
minimize('vacuum_crude')
state.gmxcalls['insert-molecules'] = {'command':'insert-molecules','required':[],'flags':[]}
state.gmxpaths['insert-molecules'] = 'gmx insert-molecules'
box_size = settings.box_size
if type(box_size) in [float,int]: box_size = ' '.join(['%s'%box_size for i in range(3)])
gmx('insert-molecules',ci='vacuum_crude-minimized',
	nmol=settings.nmol,box=box_size,log='insert-molecules',o='vacuum-melt')
component(settings.molecule_name,count=settings.nmol)
solvate(structure='vacuum-melt',gro='solvate-untrimmed')
hydration_adjust(structure='solvate-untrimmed',gro='solvate')
write_topology('solvate.top')
minimize('solvate')
copy_file('solvate-minimized.gro','system-residues.gro')
antifreeze_ratio = state.q('antifreeze_ratio',False)
if antifreeze_ratio: 
	import ipdb;ipdb.set_trace()
	add_antifreeze(structure='system-residues',gro='system-penultimate',ratio=antifreeze_ratio)
else: copy_file('system-residues.gro','system-penultimate.gro')
write_structure_by_chain(structure='system-penultimate',gro='system')
write_topology('system.top')
equilibrate()
