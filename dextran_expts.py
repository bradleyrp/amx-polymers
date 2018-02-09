{

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
	Note that this is DEPRECATED by the non-injective method in dextran_model_building.

step: tune-fine-to-coarse
mapping spec: @polymers/dextran_atomistic_to_martini.yaml
model spec: @polymers/dextran_cg_model_v1.yaml
atomistic reference:| {
	'path':'../melt-v011/s01-melt',
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
	Note that this method uses the "bigger" keyword to scale-up the melt to a size larger than the 
		atomistic source data.
	Next development is to incorporate these structures into a melt.
		see: dextran_model_building_melt

step: tune-fine-to-coarse
mapping spec: @polymers/dextran_atomistic_to_martini.yaml
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
		'input-md-short1-eq-in.mdp':[{'dt':0.001,'tau_p':1.0,'compressibility':'5e-5'}],
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
"""},

'dextran_model_building_melt':{
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
