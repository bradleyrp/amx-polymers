{

'dextran_pentamer':{
#####
####
###
##
#
'tags':['aamd'],
'script':'script-gel.py',
'params':'parameters.py',
'extensions':['codes/melts.py'],
'settings':"""

USAGE NOTES:|
	note that this was superceded first by the off-lattice, then by the tuning methods below
	DEMO: MAKE A POLYMER GEL
	SYSTEM: AGLC in CHARMM36
	needs: solvate, tune for larger systems to avoid segfaults, etc.
	2016.12.09 update:
		this script is running on dark to collect equilibrium statistics for the angle/dihedral 
			between monomers for a 5-mer
		when it's complete ryan will share the run
	plan for getting a 36-mer:
		1. finish 5-mer
		2. collect equilibrium angles
		3. decide on the "second" angle definition and collect those angles
		4. copy the code that rotates the monomers onto a linker from the 3D-lattice code
		5. write generic code for a 3D random walk with our two angles
		6. add the generic polymer 3D-walk code to the automacs script 
			(similar to the 3D lattice version) and generate a 36-mer
		7. test that the 36-mer is stable
	2017.2.28 update:
		ryan is porting this for the new version of automacs
		codes taken from melt-v012 on dark
		note there are discarded codes, viewer scripts, etc in this directory (possibly worth salvaging)
		currently porting the on-lattice copy
			note that it fails every fifth time at the vacuum step
			this is part of why we want to do the off-lattice version
		off lattice version continues below in dextran_off_lattice_36mer
		and that was later superceded by a forward-and-reverse mapping scheme because of stability problems

step: melt
polymer name: AGLC
equilibration: ""
ff_includes: []
include_adds: []
files: ['@structure-repo/dextran/aglc.gro']
sources: ['@charmm/charmm36.ff']
aglc source: @structure-repo/dextran/aglc-example.top
water buffer: 1.2
on lattice: True
review3d: True
melt settings:|{
	'a0':0.439,
	'sizer':6,
	'n_p':5,
	'volume_limit':0.99,
	'uniform':True,
	'diagonals':False,
	'review':False}

mdp specs:|{
	'group':'aamd',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep','temperature':'off'}],
		'input-md-in.mdp':[{'temperature':'bussi-other-water'}],
		}
	}

place specs:|{
	'monomer':'algc.gro',
	'repeat_unit':{'end':'O6','start':'O1','next':'C1'},
	'linkage_delete_atoms':[['HO6'],['HO1','O1']],
	'atom_name_changes':[
		{'from':'OC311','from_charge':-0.65,'to_charge':-0.36,
		'to':'OC301','atom':'O6','which':'previous'},
		{'from':'CC3162','from_charge':0.34,'to_charge':0.29,
		'to':'CC3162','atom':'C1','which':'next'},
		{'from':'CC321','from_charge':0.050,'to_charge':0.00,
		'to':'CC321','atom':'C6','which':'previous'}]}

"""},

'dextran_off_lattice_36mer':{
#####
####
###
##
#
'tags':['aamd'],
'script':'script-gel.py',
'params':'parameters.py',
'extensions':['codes/melts.py'],
'settings':"""

USAGE NOTES:|
	copied from dextran_dev_pentamer
	note that this was superceded by the tuning method below

step: melt
polymer name: AGLC
equilibration: ""
force field: charmm36
files: ['@structure-repo/dextran/aglc.gro']
sources: ['@charmm/charmm36.ff']
aglc source: @structure-repo/dextran/aglc-example.top
water buffer: 1.2
on lattice: False
review3d: False
melt settings:|{
	'a0':0.439,
	'sizer':6,
	'n_p':36,
	'volume_limit':0.99,
	'uniform':True,
	'diagonals':False,
	'review':False,
	'angle':-145.0,
	'torsion':100.0}

mdp specs:|{
	'group':'aamd',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep','temperature':'off'}],
		'input-md-in.mdp':[{'temperature':'bussi-other-water'}],},}

place specs:|{
	'monomer':'algc.gro',
	'repeat_unit':{'end':'O6','start':'O1','next':'C1'},
	'linkage_delete_atoms':[['HO6'],['HO1','O1']],
	'atom_name_changes':[
		{'from':'OC311','from_charge':-0.65,'to_charge':-0.36,
		'to':'OC301','atom':'O6','which':'previous'},
		{'from':'CC3162','from_charge':0.34,'to_charge':0.29,
		'to':'CC3162','atom':'C1','which':'next'},
		{'from':'CC321','from_charge':0.050,'to_charge':0.00,
		'to':'CC321','atom':'C6','which':'previous'}]}

"""},

'dextran_coarse_basic':{
#####
####
###
##
#
'tags':['cgmd','tested_2017.09.14'],
'script':'script-cg-gel.py',
'params':'parameters.py',
'extensions':['codes/melts.py'],
'settings':"""

USAGE NOTES:|
	copied from dextran_dev_off_lattice_36mer
	goal is to simulation a crude version of martini poly-glucose chain (maybe without bimodal torsions)
		as a starting point for backmapping the atomistic model of dextran
	note that this simulation is not always stable

step: melt
polymer name: AGLC
equilibration: ['short']
force field: martini-sources
aglc source: None
water buffer: 8
on lattice: False
review3d: False
sol: W
solvent: martini-water
melt settings:|{
	'n_p':36,
	'a0':0.356,
	'angle':90.0,
	'torsion':142.0,}

files: ['@martini/library-general-structs/martini-water.gro']
sources: ['inputs/martini/martini-sources.ff']

mdp specs:|{
	'group':'cgmd-polymers',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep'}],
		'input-md-short-eq-in.mdp':[{'dt':0.001,'tau_p':1.0,'compressibility':'5e-5'}],
		'input-md-in.mdp':[{'dt':0.02}],},}

place specs:|{
	'monomer':'algc.gro',
	'repeat_unit':{'end':'O6','start':'O1','next':'C1'},
	'linkage_delete_atoms':[['HO6'],['HO1','O1']],
	'atom_name_changes':[
		{'from':'OC311','from_charge':-0.65,'to_charge':-0.36,
		'to':'OC301','atom':'O6','which':'previous'},
		{'from':'CC3162','from_charge':0.34,'to_charge':0.29,
		'to':'CC3162','atom':'C1','which':'next'},
		{'from':'CC321','from_charge':0.050,'to_charge':0.00,
		'to':'CC321','atom':'C6','which':'previous'}]}

"""},

'dextran_coarse_tuner_basic':{
#####
####
###
##
#
'tags':['cgmd','dev'],
'script':'script-tune-cg.py',
'params':'parameters.py',
'extensions':['codes/melts.py','codes/melts_tuner.py'],
'settings':"""
step: tune
method style: basic

USAGE NOTES:|
	temporary placeholder for tuning a new MARTINI model for dextran based on the pentamer in CHARMM

"""},

'dextran_coarse_tuner':{
#####
####
###
##
#
'tags':['cgmd','dev'],
'script':'script-tune-cg.py',
'params':'parameters.py',
'extensions':['codes/melts.py','codes/melts_tuner.py'],
'settings':"""
step: tune
method style: advanced

USAGE NOTES:|
	temporary placeholder for tuning a new MARTINI model for dextran based on the pentamer in CHARMM

"""},

'dextran_scale_switching':{
#####
####
###
##
#
'tags':['aamd_cgmd','dev'],
'metarun':[
{'step':'gel','do':'dextran_coarse_basic','settings':"""

USAGE NOTES:|
	created on 2017.07.06.1400
	summary of work to date:
		dextran_coarse_tuner (formerly dextran_martini_dev_tuner)
			reads in the pentamer in v011
			has a section at the end that backmaps
				from one v011 frame
				to another v011 frame that has been coarse-grained
			and hence represents the simplest possible backmapping we could do
		dextran_coarse_basic (formerly dextran_martini_dev)
			make a larger N-mer of dextran Martini
			with reasonable guesses about the various parameters
	purpose of this experiment:
		1. make the gel with best-guess parameters
		2. take a frame from that gel and backmap the atomistic coordinates onto them
		3. test for stability in aamd? if unstable, then try iterating.

water buffer: 2
melt settings:|{
	'n_p':8,
	'a0':0.356,
	'angle':90.0,
	'torsion':142.0,}

"""},]},

'dextran_backmapping':{
#####
####
###
##
#
'tags':['aamd_cgmd','dev'],
'metarun':[
{'step':'gel','do':'dextran_coarse_tuner_basic'},
{'step':'gel','do':'dextran_coarse_tuner'},]},

### SECOND ROUND OF DEVELOPMENT, stored here for posterity on 2018.02.08

'dextran_coarse_tuner':{
#####
####
###
##
#
'tags':['cgmd','dev'],
'script':'script-tune-cg.py',
'params':'parameters.py',
'extensions':['codes/melts.py'],
'settings':"""
# this experiment is a placeholder
step: tune
method style: backward
mapping spec: @polymers/dextran_atomistic_to_martini.yaml
"""},

'dextran_model_building_transitional':{
#####
####
###
##
#
'tags':['aamd_cgmd','dev'],
'metarun':[
{'step':'gel','do':'dextran_coarse_tuner','settings':"""
method style: forward
step: tune-forwards
USE NOTES:|
	You have to prepare the trajectory correctly:
	gmx trjconv -s md.part0005.tpr -f md.parts.pbcmol.xtc -o md.parts.pbcmol.centered.xtc -fit progressive
atomistic reference:| {
	'path':'/home/rpb/omicron/dataset-project-polymers/melt-v011/s01-melt',
	'gro':'system.gro','xtc':'md.parts.pbcmol.centered.xtc','n_monomers':5,
	'molecule_map':'one molecule many residues injective','selection':'resname AGLC'}
"""},
{'step':'gel','do':'dextran_coarse_tuner','settings':"""
step: tune-backwards
method style: backward
coarse reference:| {
	'path':'/home/rpb/omicron/dataset-project-polymers/melt-v027-dextran-scale-switching-coarse-8mer/s01-melt',
	'gro':'md.part0001.gro','n_monomers':9,}
atomistic reference:| {
	'path':'/home/rpb/omicron/dataset-project-polymers/melt-v011/s01-melt',
	'gro':'system.gro','xtc':'md.parts.pbcmol.xtc','n_monomers':5,}
"""},]},

'dextran_coarse_basic':{
#####
####
###
##
#
'tags':['cgmd','tested_2017.09.14'],
'script':'script-cg-gel.py',
'params':'parameters.py',
'extensions':['codes/melts.py'],
'settings':"""

USAGE NOTES:|
	copied from dextran_dev_off_lattice_36mer
	goal is to simulation a crude version of martini poly-glucose chain (maybe without bimodal torsions)
		as a starting point for backmapping the atomistic model of dextran
	note that this simulation is not always stable

step: melt
polymer name: AGLC
equilibration: ['short']
force field: martini-sources
aglc source: None
water buffer: 8
on lattice: False
review3d: False
sol: W
solvent: martini-water
melt settings:|{
	'n_p':36,
	'a0':0.356,
	'angle':90.0,
	'torsion':142.0,}

files: ['@martini/library-general-structs/martini-water.gro']
sources: ['inputs/martini/martini-sources.ff']

mdp specs:|{
	'group':'cgmd-polymers',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep'}],
		'input-md-short-eq-in.mdp':[{'dt':0.001,'tau_p':1.0,'compressibility':'5e-5'}],
		'input-md-in.mdp':[{'dt':0.02}],},}

place specs:|{
	'monomer':'algc.gro',
	'repeat_unit':{'end':'O6','start':'O1','next':'C1'},
	'linkage_delete_atoms':[['HO6'],['HO1','O1']],
	'atom_name_changes':[
		{'from':'OC311','from_charge':-0.65,'to_charge':-0.36,
		'to':'OC301','atom':'O6','which':'previous'},
		{'from':'CC3162','from_charge':0.34,'to_charge':0.29,
		'to':'CC3162','atom':'C1','which':'next'},
		{'from':'CC321','from_charge':0.050,'to_charge':0.00,
		'to':'CC321','atom':'C6','which':'previous'}]}

"""},

'dextran_model_building_injective':{
#####
####
###
##
#
'tags':['cgmd','dev'],
'script':'script-tune-cg.py',
'params':'parameters.py',
'extensions':['codes/melts.py'],
'settings':"""

USE NOTES:|
	Master experiment which analyzes an atomistic trajectory, builds and runs a coarse-model,
		and then attempts a backmapping.
	When this is complete, we will add it to a metarun that creates that atomistic data.
	This was spun-off from dextran_model_building_transitional above.
	You have to prepare the trajectory correctly:
		gmx trjconv -s md.part0005.tpr -f md.parts.pbcmol.xtc -o md.parts.pbcmol.centered.xtc -fit progressive

step: tune-fine-to-coarse
mapping spec: @polymers/dextran_atomistic_to_martini.yaml
model spec: @polymers/dextran_cg_model_v1.yaml
atomistic reference:| {
	'path':'/home/rpb/omicron/dataset-project-polymers/melt-v011/s01-melt',
	'gro':'system.gro','xtc':'md.parts.pbcmol.centered.xtc','n_monomers':5,
	'molecule_map':'one molecule many residues injective','selection':'resname AGLC'}
#! either use injective above and n_p 5 below or bigger above and something bigger for n_p
# additions for CG model
files: ['@martini/library-general-structs/martini-water.gro']
sources: ['inputs/martini/martini-sources.ff']
mdp specs:|{
	'group':'cgmd-polymers',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep'}],
		'input-md-short-eq-in.mdp':[{'dt':0.001,'tau_p':1.0,'compressibility':'5e-5'}],
		'input-md-in.mdp':[{'dt':0.02}],},}
place specs:|{
	'monomer':'algc.gro',
	'repeat_unit':{'end':'O6','start':'O1','next':'C1'},
	'linkage_delete_atoms':[['HO6'],['HO1','O1']],
	'atom_name_changes':[
		{'from':'OC311','from_charge':-0.65,'to_charge':-0.36,
		'to':'OC301','atom':'O6','which':'previous'},
		{'from':'CC3162','from_charge':0.34,'to_charge':0.29,
		'to':'CC3162','atom':'C1','which':'next'},
		{'from':'CC321','from_charge':0.050,'to_charge':0.00,
		'to':'CC321','atom':'C6','which':'previous'}]}
polymer name: AGLC
equilibration: ['short']
force field: martini-sources
aglc source: None
water buffer: 8
on lattice: False
review3d: False
sol: W
solvent: martini-water
#! retaining angle/torsion for building on the lattice
#! note that n_p must match the incoming source if you are using "one molecule many residues injective"
melt settings:|{
	'n_p':5,
	'a0':0.356,
	'angle':90.0,
	'torsion':142.0,}

"""},

'dextran_model_building':{
#####
####
###
##
#
'tags':['cgmd','dev'],
'script':'script-tune-cg.py',
'params':'parameters.py',
'extensions':['codes/melts.py'],
'settings':"""

USE NOTES:|
	Master experiment which analyzes an atomistic trajectory, builds and runs a coarse-model,
		and then attempts a backmapping.
	When this is complete, we will add it to a metarun that creates that atomistic data.
	This was spun-off from dextran_model_building_transitional above.
	You have to prepare the trajectory correctly:
		gmx trjconv -s md.part0005.tpr -f md.parts.pbcmol.xtc -o md.parts.pbcmol.centered.xtc -fit progressive

step: tune-fine-to-coarse
mapping spec: @polymers/dextran_atomistic_to_martini.yaml
model spec: @polymers/dextran_cg_model_v1.yaml
atomistic reference:| {
	'path':'/home/rpb/omicron/dataset-project-polymers/melt-v011/s01-melt',
	'gro':'system.gro','xtc':'md.parts.pbcmol.centered.xtc','n_monomers':5,
	'molecule_map':'one molecule many residues bigger','selection':'resname AGLC'}
#! either use injective above and n_p 5 below or bigger above and something bigger for n_p
# additions for CG model
files: ['@martini/library-general-structs/martini-water.gro']
sources: ['inputs/martini/martini-sources.ff']
mdp specs:|{
	'group':'cgmd-polymers',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep'}],
		'input-md-short1-eq-in.mdp':[{'dt':0.001,'Pcoupl':'no'}],
		'input-md-short2-eq-in.mdp':[{'dt':0.001,'tau_p':1.0,'compressibility':'5e-5'}],
		'input-md-short3-eq-in.mdp':[{'dt':0.005,'tau_p':1.0,'compressibility':'5e-5'}],
		'input-md-in.mdp':[{'dt':0.02}],},}
place specs:|{
	'monomer':'algc.gro',
	'repeat_unit':{'end':'O6','start':'O1','next':'C1'},
	'linkage_delete_atoms':[['HO6'],['HO1','O1']],
	'atom_name_changes':[
		{'from':'OC311','from_charge':-0.65,'to_charge':-0.36,
		'to':'OC301','atom':'O6','which':'previous'},
		{'from':'CC3162','from_charge':0.34,'to_charge':0.29,
		'to':'CC3162','atom':'C1','which':'next'},
		{'from':'CC321','from_charge':0.050,'to_charge':0.00,
		'to':'CC321','atom':'C6','which':'previous'}]}
polymer name: AGLC
equilibration: ['short1','short2','short3']
force field: martini-sources
aglc source: None
water buffer: 8
on lattice: True
review3d: False
sol: W
solvent: martini-water
#! retaining angle/torsion for building on the lattice
#! note that n_p must match the incoming source if you are using "one molecule many residues injective"
melt settings:|{
	'n_p':30,
	'a0':0.356,
	'angle':90.0,
	'torsion':142.0,
	'volume_limit':0.5}
#! the following settings override the ones above
lattice melt settings:|{
	'n_p':30,
	'volume_limit':0.05,
	'a0':0.35,
	'sizer':20}
"""},

}