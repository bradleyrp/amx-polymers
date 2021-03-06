{

'maltoheptaose':{
#####
####
###
##
#
'tags':['cgmd','test'],
'script':'script-maltoheptaose.py',
'params':'parameters.py',
'extensions':['melts.py','melts_simple.py'],
'settings':"""

USAGE NOTES:|
	maltoheptaose melt using standard parameters
	this experiment was designed to make a stock martini oligo carbohydrate for comparison
	we demoed the oligoglucose_coarse_simple which used the random walk method to make one oligomer
	this method seeks a melt, using the stock martini parameters instead of the first guess

step:               melt              # folder name for this step
n_p:                7                 # number of monomers in a polymer
molecule name:      Maltoheptaose     # molecule name in the final topology
residue name:       AMYL              # residue name in the molecule
a0:                 0.35              # backbone spacer for initial guess
water ratio:        10                # water beads per polymer beads

melt build method:  insert-molecules  # use the standard gromacs routine for making the melt
beads per monomer:  3                 # beads in each monomer
backbone position:  3                 # which bead is the backbone
nmol:               100               # molecules for `gmx insert-molecules`
box size:           12.               # box size for `gmx insert-molecules`
antifreeze ratio:   0.0               # whether to add antifreeze

# a string holding a lambda function for naming beads
#! consider modifying this to read the ITP directly?
atom namer: "lambda x: 'B%d'%([2,3,1,4,6,5][x%6])"
#! hack below
rename_detected_composition: {'AMYL':'Maltoheptaose'}
composition_adjust: "def composition_adjust(composition):\\n\\tcomposition_adjusted = []\\n\\tfor i,j in composition:\\n\\t\\tcomposition_adjusted.append([i,int(j/7.) if i=='AMYL' else j])\\n\\treturn composition_adjusted"
maxwarn: 4

sol: W
solvent: martini-water
force field: martini-sources
sources: ['@martini/martini-sources.ff']
files:| [
	'@martini/martini-extensions/martini_v2.0_sugars.itp',
	'@martini/library-general-structs/martini-water.gro']

equilibration: ['short1','short2']
mdp specs:|{
	'group':'cgmd-polymers',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep'}],
		'input-md-short1-eq-in.mdp':[{'dt':0.001,'Pcoupl':'no'}],
		'input-md-short2-eq-in.mdp':[{'dt':0.001,'ref_p':'1.0 1.0 1.0 0.0 0.0 0.0','tau_p':'1.0','Pcoupltype':'anisotropic','compressibility':'5e-5 5e-5 5e-5 0.0 0.0 0.0'}],
		'input-md-short3-eq-in.mdp':[{'dt':0.005,'ref_p':'1.0 1.0 1.0 0.0 0.0 0.0','tau_p':'1.0','Pcoupltype':'anisotropic','compressibility':'5e-5 5e-5 5e-5 0.0 0.0 0.0'}],
		'input-md-in.mdp':[{'dt':0.02}],},}

"""},

'maltoheptaose_big':{
#####
####
###
##
#
'tags':['cgmd','tested_2018.04.27'],
'metarun':[
{'step':'small','do':'maltoheptaose','settings':"""
nmol:               10  # molecules for `gmx insert-molecules`
box size:           5.  # box size for `gmx insert-molecules`
#! hack below
rename_detected_composition: {'AMYL':'Maltoheptaose'}
composition_adjust: "def composition_adjust(composition):\\n\\tcomposition_adjusted = []\\n\\tfor i,j in composition:\\n\\t\\tcomposition_adjusted.append([i,int(j/7.) if i=='AMYL' else j])\\n\\treturn composition_adjusted"
water ratio: 10.
antifreeze ratio: 0.
"""},
{'step':'large','do':'multiply_general','settings':"""
step: large
requires: multiply
equilibration: ['short1','short2','short3']
#! repeated for restarts and write_structure_by_chain
n_p:                7                 # number of monomers in a polymer
molecule name:      Maltoheptaose     # molecule name in the final topology
residue name:       AMYL              # residue name in the molecule
beads per monomer:  3                 # beads in each monomer
#! hack below
composition_adjust: "def composition_adjust(composition):\\n\\tcomposition_adjusted = []\\n\\tfor i,j in composition:\\n\\t\\tcomposition_adjusted.append([i,int(j/7.) if i=='AMYL' else j])\\n\\treturn composition_adjusted"
rename_detected_composition: {'AMYL':'Maltoheptaose'}
maxwarn: 4
minimize: True
proceed: True
genconf gap: 0.3
nx: 2
ny: 2
nz: 2
"""}]},

'glycans':{
#####
####
###
##
#
'tags':['cgmd','tested_2018.05.01'],
'metarun':[
{'do':'maltoheptaose','settings':"""

molecule name: GLUC
residue name: GLUC
n_p: 1
water ratio: 10
nmol: 40 
box size: 5.
antifreeze ratio: 0
atom namer: "lambda x: 'B%d'%([1,2,3][x])"

"""},
{'step':'large-melt','do':'multiply_general','settings':"""
step: large
requires: multiply
equilibration: ['short1','short2','short3']
rename_detected_composition: {'GLUC':'GLUC'}
#! repeated for restarts and write_structure_by_chain
molecule name: GLUC
residue name: GLUC
n_p: 1
beads per monomer: 3
#! end hack
maxwarn: 4
minimize: True
proceed: True
genconf gap: 0.3
nx: 2
ny: 2
nz: 2
"""}]},

'amylose':{
#####
####
###
##
#
'tags':['cgmd','test'],
'script':'script-amylose.py',
'params':'parameters.py',
'extensions':['melts.py','melts_simple.py','amylose.py'],
'settings':"""

USAGE NOTES:|
	maltoheptaose melt using standard parameters
	this experiment was designed to make a stock martini oligo carbohydrate for comparison
	we demoed the oligoglucose_coarse_simple which used the random walk method to make one oligomer
	this method seeks a melt, using the stock martini parameters instead of the first guess

step:               melt              # folder name for this step
n_p:                41                # number of monomers in a polymer
molecule name:      AmylBig           # molecule name in the final topology
residue name:       AMYL              # residue name in the molecule
a0:                 0.35              # backbone spacer for initial guess
water ratio:        5                 # water beads per polymer beads

melt build method:  insert-molecules  # use the standard gromacs routine for making the melt
beads per monomer:  3                 # beads in each monomer
backbone position:  3                 # which bead is the backbone
nmol:               10                # molecules for `gmx insert-molecules`
box size:           8.                # box size for `gmx insert-molecules`
antifreeze ratio:   0.0               # whether to add antifreeze
do weird solvate:   False             # weird solvate attempts to solve memory/time problems on big systems

# source topology
heptaose source: @martini/martini-extensions/martini_v2.0_sugars.itp
# a string holding a lambda function for naming beads
#! consider modifying this to read the ITP directly?
atom namer: "lambda x: 'B%d'%([2,3,1,4,6,5][x%6])"
#! hack below
rename_detected_composition: {'AMYL':'AmylBig'}
composition_adjust: "def composition_adjust(composition):\\n\\tcomposition_adjusted = []\\n\\tfor i,j in composition:\\n\\t\\tcomposition_adjusted.append([i,int(j/7.) if i=='AMYL' else j])\\n\\treturn composition_adjusted"
maxwarn: 4

water gap: 5

sol: W
solvent: martini-water
force field: martini-sources
sources: ['@martini/martini-sources.ff']
files:| [
	'@martini/martini-extensions/martini_v2.0_sugars.itp',
	'@martini/library-general-structs/martini-water.gro']

equilibration: ['short1','short2','short3']
mdp specs:|{
	'group':'cgmd-polymers',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep'}],
		'input-md-short1-eq-in.mdp':[{'dt':0.0001,'Pcoupl':'no'}],
		'input-md-short2-eq-in.mdp':[{'dt':0.001}],
		'input-md-short3-eq-in.mdp':[{'dt':0.005,'nsteps':2000000}],
		'input-md-in.mdp':[{'dt':0.02}],},}
# removed from shor2 short3: ,'ref_p':'1.0 1.0 1.0 0.0 0.0 0.0','tau_p':'1.0','Pcoupltype':'anisotropic','compressibility':'5e-5 5e-5 5e-5 0.0 0.0 0.0'

"""},

'amylose_big':{
#####
####
###
##
#
'tags':['cgmd','tested_2018.04.27'],
'metarun':[
{'step':'small','do':'amylose','settings':"""

step:               melt              # folder name for this step
n_p:                41                # number of monomers in a polymer
molecule name:      AmylBig           # molecule name in the final topology
residue name:       AMYL              # residue name in the molecule
a0:                 0.35              # backbone spacer for initial guess
water ratio:        5                 # water beads per polymer beads

melt build method:  insert-molecules  # use the standard gromacs routine for making the melt
beads per monomer:  3                 # beads in each monomer
backbone position:  3                 # which bead is the backbone
nmol:               10                # molecules for `gmx insert-molecules`
box size:           10.               # box size for `gmx insert-molecules`
antifreeze ratio:   0.0               # whether to add antifreeze
do weird solvate:   False             # weird solvate attempts to solve memory/time problems on big systems

# source topology
heptaose source: @martini/martini-extensions/martini_v2.0_sugars.itp
# a string holding a lambda function for naming beads
#! consider modifying this to read the ITP directly?
atom namer: "lambda x: 'B%d'%([2,3,1,4,6,5][x%6])"
#! hack below
rename_detected_composition: {'AMYL':'AmylBig'}
composition_adjust: "def composition_adjust(composition):\\n\\tcomposition_adjusted = []\\n\\tfor i,j in composition:\\n\\t\\tcomposition_adjusted.append([i,int(j/7.) if i=='AMYL' else j])\\n\\treturn composition_adjusted"
maxwarn: 4
water gap: 5
sol: W
solvent: martini-water
force field: martini-sources
sources: ['@martini/martini-sources.ff']
files:| [
	'@martini/martini-extensions/martini_v2.0_sugars.itp',
	'@martini/library-general-structs/martini-water.gro']

equilibration: ['short1','short2','short3']
mdp specs:|{
	'group':'cgmd-polymers',
	'mdps':{
		'input-em-steep-in.mdp':[{'integrator':'steep'}],
		'input-md-short1-eq-in.mdp':[{'dt':0.0001,'Pcoupl':'no'}],
		'input-md-short2-eq-in.mdp':[{'dt':0.001}],
		'input-md-short3-eq-in.mdp':[{'dt':0.005,'nsteps':2000000}],
		'input-md-in.mdp':[{'dt':0.02}],},}
# removed from shor2 short3: ,'ref_p':'1.0 1.0 1.0 0.0 0.0 0.0','tau_p':'1.0','Pcoupltype':'anisotropic','compressibility':'5e-5 5e-5 5e-5 0.0 0.0 0.0'

"""},
{'step':'large','do':'multiply_general','settings':"""
step: large
requires: multiply
equilibration: ['short1','short2','short3']
#! repeated for restarts and write_structure_by_chain
n_p:                7                 # number of monomers in a polymer
molecule name:      Maltoheptaose     # molecule name in the final topology
residue name:       AMYL              # residue name in the molecule
beads per monomer:  3                 # beads in each monomer
#! hack below
composition_adjust: "def composition_adjust(composition):\\n\\tcomposition_adjusted = []\\n\\tfor i,j in composition:\\n\\t\\tcomposition_adjusted.append([i,int(j/7.) if i=='AMYL' else j])\\n\\treturn composition_adjusted"
rename_detected_composition: {'AMYL':'Maltoheptaose'}
maxwarn: 4
minimize: True
proceed: True
genconf gap: 0.3
nx: 3
ny: 3
nz: 3
"""}]},

}
