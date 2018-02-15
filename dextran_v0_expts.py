{


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

'1dextran_coarse_tuner':{
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

'2dextran_coarse_tuner':{
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

}