#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
write_mdp()
component(settings.molecule_name,count=1)
state.itp = [settings.itp_file]
state.gmxcalls['insert-molecules'] = {'command':'insert-molecules','required':[],'flags':[]}
state.gmxpaths['insert-molecules'] = 'gmx insert-molecules'
box_size = settings.box_size
if type(box_size) in [float,int]: box_size = ' '.join(['%s'%box_size for i in range(3)])
gmx('insert-molecules',ci=settings.structure,
	nmol=settings.nmol,box=box_size,log='insert-molecules',o='vacuum')
component(settings.molecule_name,count=settings.nmol)
solvate(structure='vacuum',gro='solvate')
write_topology('solvate.top')
minimize('solvate')
counterions(
	structure='solvate-minimized',
	top='solvate',
	ff_includes='ions')
minimize('counterions')
write_structure_pdb(
	pdb='start-structure.pdb',
	structure='counterions')
write_top('system.top')
write_continue_script()
equilibrate()
