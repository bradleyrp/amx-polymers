{

'dextran_dev_pentamer':{
#####
####
###
##
#
'tags':['aamd','tag_dev'],
'script':'script-gel.py',
'params':'parameters.py',
'extensions':['codes/melts.py'],
'settings':"""

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

USAGE NOTES:|
	DEMO: MAKE A POLYMER GEL
	SYSTEM: AGLC in CHARMM36
	needs: solvate, tune for larger systems to avoid segfaults, etc.
	2016.12.09 update:
		this script is running on dark to collect equilibrium statistics for the angle/dihedral 
			between monomers for a 5-mer
		when it's complete ryan will upload to green
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
"""},

'dextran_dev_off_lattice_36mer':{
#####
####
###
##
#
'tags':['aamd','tag_dev'],
'script':'script-gel.py',
'params':'parameters.py',
'extensions':['codes/melts.py'],
'settings':"""

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

USAGE NOTES:|
	copied from dextran_dev_pentamer
"""},

'dextran_martini_dev':{
#####
####
###
##
#
'tags':['cgmd','tag_dev','tested_2017.09.14'],
'script':'script-cg-gel.py',
'params':'parameters.py',
'extensions':['codes/melts.py'],
'settings':"""

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
	'torsion':142.0,
	}

files: ['@martini/library-general-structs/martini-water.gro']
sources:| [
    'inputs/martini/martini-sources.ff',
    ]

mdp specs:|{
	'group':'cgmd-polymers',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep'}],
		'input-md-short-eq-in.mdp':[{'dt':0.001,'tau_p':1.0,'compressibility':'5e-5'}],
		'input-md-in.mdp':[],
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

USAGE NOTES:|
	copied from dextran_dev_off_lattice_36mer
	goal is to simulation a crude version of martini poly-glucose chain (maybe without bimodal torsions)
		as a starting point for backmapping the atomistic model of dextran
"""},

'dextran_martini_dev_tuner':{
#####
####
###
##
#
'tags':['cgmd','tag_dev'],
'script':'script-tune-cg.py',
'params':'parameters.py',
'extensions':['codes/melts.py','codes/melts_tuner.py'],
'settings':"""

step: tune

USAGE NOTES:|
	THIS IS A TEMPORARY PLACEHOLDER TO TUNE A NEW MARTINI MODEL FOR DEXTRAN FROM THE PENTAMER IN CHARMM
"""},

'dextran_scale_switching':{
#####
####
###
##
#
'tags':['tag_dev','aamd_cgmd'],
'metarun':[
{'step':'gel','do':'dextran_martini_dev','settings':"""

water buffer: 2
melt settings:|{
	'n_p':8,
	'a0':0.356,
	'angle':90.0,
	'torsion':142.0,
	}

USAGE NOTES:|
	created on 2017.07.06.1400
	summary of work to date:
		dextran_martini_dev_tuner
			reads in the pentamer in v011
			has a section at the end that backmaps
				from one v011 frame
				to another v011 frame that has been coarse-grained
			and hence represents the simplest possible backmapping we could do
		dextran_martini_dev
			make a larger N-mer of dextran Martini
			with reasonable guesses about the various parameters
	purpose of this experiment:
		1. make the gel with best-guess parameters
		2. take a frame from that gel and backmap the atomistic coordinates onto them
		3. test for stability in aamd? if unstable, then try iterating.
"""},]}

}