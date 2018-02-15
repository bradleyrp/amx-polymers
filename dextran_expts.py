{

### SIMPLE METHODS

'glucose_16_atomistic_pentamer':{
#####
####
###
##
#
'tags':['aamd','tested_2018.02.13'],
'script':'script-polymer-simple.py',
'params':'parameters.py',
'extensions':['melts_simple.py'],
'settings':"""

USAGE NOTES:|
	this procedure makes a glucose pentamer with 1,6 linkages using a starting monomer structure
	uses CHARMM36 which generated the model and starting topology
	if you extend this beyond an octamer, there are stability problems
	this simulation serves as source material for other experiments
	stability problems when making large atomistic systems motivated our implementation of a CG model

step: melt
polymer name: AGLC
equilibration: ""
ff_includes: []
include_adds: []
files: ['@polymers/starts/aglc.gro']
sources: ['@charmm/charmm36.ff']
aglc source: @polymers/starts/aglc-example.top
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

'oligoglucose_coarse_simple':{
#####
####
###
##
#
'tags':['cgmd','tested_2017.02.11'],
'script':'script-polymer-simple-coarse.py',
'params':'parameters.py',
'extensions':['melts_simple.py'],
'settings':"""

USAGE NOTES:|
	a simple version of a coarse-grained oligo-glucose model
	creates the oligomer off-lattice
	note that subsequent changes to parameters.py may cause this to freeze

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

### ADVANCED METHODS

'tune_polymer_pentamer':{
#####
####
###
##
#
'tags':['cgmd','tested_2018.02.14'],
'script':'script-multiscale.py',
'params':'parameters.py',
'extensions':['melts.py','melts_simple.py','plotter_omni_panels.py'],
'settings':"""

USE NOTES:|
	This experiment was the first to use the MultiScaleModel class

	Master experiment which analyzes an atomistic trajectory, then builds and runs a coarse-model.
	When this is complete, we will add it to a metarun that creates that atomistic data.
	You have to prepare the trajectory correctly:
		gmx trjconv -s md.part0005.tpr -f md.parts.pbcmol.xtc \
			-o md.parts.pbcmol.centered.xtc -fit progressive
	Note that this is SUPERCEDED by the non-injective method in dextran_model_building.
	DEVELOPMENT NOTE: currently crashes on part0001 due to pressure problems.
		this is possibly a result of the small system size. see dextran_model_building_single_large below

step: tune-fine-to-coarse
mapping spec: @polymers/dextran_atomistic_to_martini_v1.yaml
model spec: @polymers/dextran_cg_model_v1.yaml
atomistic reference:| {
	'path':'../melt-v011/s01-melt',
	'gro':'system.gro','xtc':'md.parts.pbcmol.centered.xtc','n_monomers':5,
	'molecule_map':'one molecule many residues injective','selection':'resname AGLC'}
files: ['@martini/library-general-structs/martini-water.gro']
sources: ['inputs/martini/martini-sources.ff']
maxwarn: 2
mdp specs:|{
	'group':'cgmd-polymers',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep'}],
		'input-md-short1-eq-in.mdp':[{'dt':0.001}],
		'input-md-short2-eq-in.mdp':[{'dt':0.01}],
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
equilibration: ['short1','short2']
force field: martini-sources
aglc source: None
water buffer: 8
on lattice: False
review3d: False
sol: W
solvent: martini-water
#! note that n_p must match the source data for the injective method
#! note that the angle and torsion set the initial guess for the structure
melt settings: {'n_p':5,'a0':0.356,'angle':90.0,'torsion':142.0,}
"""},

'tune_polymer_large_single':{
#####
####
###
##
#
'tags':['cgmd','tested_2018.02.14'],
'script':'script-multiscale.py',
'params':'parameters.py',
'extensions':['melts.py','melts_simple.py','plotter_omni_panels.py'],
'settings':"""

USE NOTES:|
	Master experiment which analyzes an atomistic trajectory then builds a single, much bigger coarse polymer.
	This extends detran_model_building_injective by generalizing the distributions from the atomistic
		source data so that they can be mapped onto a much larger polymer.
	This method will be SUPERCEDED by the dextran_model_building_melt which uses multiple large polymers.
	Tested provisionally, with dihedrals disabled in melts.py. Stability problems need to be overcome.

step: tune-fine-to-coarse
mapping spec: @polymers/dextran_atomistic_to_martini_v1.yaml
model spec: @polymers/dextran_cg_model_v1.yaml
atomistic reference:| {
	'path':'/home/rpb/omicron/dataset-project-polymers/melt-v011/s01-melt',
	'gro':'system.gro','xtc':'md.parts.pbcmol.centered.xtc','n_monomers':5,
	'molecule_map':'one molecule many residues bigger','selection':'resname AGLC'}
files: ['@martini/library-general-structs/martini-water.gro']
sources: ['inputs/martini/martini-sources.ff']
maxwarn: 2
mdp specs:|{
	'group':'cgmd-polymers',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep'}],
		'input-md-short1-eq-in.mdp':[{'dt':0.001}],
		'input-md-short2-eq-in.mdp':[{'dt':0.01}],
		'input-md-in.mdp':[{'dt':0.02,'couple':'couple','Pcoupl':'no'}],},}
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
equilibration: ['short1','short2']
force field: martini-sources
aglc source: None
water buffer: 8
on lattice: False
review3d: False
sol: W
solvent: martini-water
# set the desired polymer size below
melt settings: {'n_p':30,'a0':0.356,'angle':90.0,'torsion':142.0}
"""},

'dextran_model_building_melt':{
#####
####
###
##
#
'tags':['cgmd','dev'],
'script':'script-multiscale.py',
'params':'parameters.py',
'extensions':['melts.py','melts_simple.py','plotter_omni_panels.py'],
'settings':"""

USE NOTES:|
	Builds on dextran_model_building_single_large by making a large melt.
	Tested provisionally, with no ...

step: tune-fine-to-coarse
mapping spec: @polymers/dextran_atomistic_to_martini_v1.yaml
model spec: @polymers/dextran_cg_model_v1.yaml
atomistic reference:| {
	'path':'../melt-v011/s01-melt',
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
# set the polymer size below
melt settings: {'n_p':30}
# settings for the starting melt structure
lattice melt settings: {'n_p':30,'volume_limit':0.05,'a0':0.35,'sizer':20}
"""},

}
