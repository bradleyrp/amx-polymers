{

'dextran_dev':{
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

step: melt
polymer name: AGLC
equilibration: ""
ff_includes: []
include_adds: []
files: ['@structure-repo/dextran/aglc.gro']
sources: ['@charmm/charmm36.ff']
aglc source: @structure-repo/dextran/aglc-example.top
water buffer: 1.2
on lattice: False
review 3d: True
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

}