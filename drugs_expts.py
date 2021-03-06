{
'doxorubicin':{
#####
####
###
##
#
'tags':['aamd','test'],
'script':'script-doxorubicin.py',
'params':'@polymers/parameters.py',
'extensions':[],
'settings':"""
step: water
molecule_name: lig
structure: doxorubicin.pdb
itp_file: doxorubicin.itp
force field: charmm36
files: ['@polymers/doxorubicin/doxorubicin.itp','@polymers/doxorubicin/doxorubicin.pdb']
sources: ['@polymers/charmm36.ff']
nmol: 1
box_size: 4
ionic strength: 0.150               # desired molar ionic strength
cation: NA                          # name of the cation for neutralizing the system
anion: CL                           # name of the anion for neutralizing the system
mdp_specs:| {
	'group':'aamd',
	'mdps':{
		'input-em-steep-in.mdp':['minimize'],
		'input-em-cg-in.mdp':['minimize',{'integrator':'cg'}],
		'input-md-in.mdp':['npt-ligand',{'nsteps':100000,}],
		},
	}
"""},
'doxorubicin_martini':{
#####
####
###
##
#
'tags':['aamd','test'],
'script':'script-coarse-small-molecule.py',
'params':'@polymers/parameters.py',
'extensions':['coarse_grain_simple.py'],
'settings':"""
step: coarse-grain
source_xtc: /home/rpb/omicron/dataset-project-polymers/melt-v048-dox/a02-aside/source.xtc
source_gro: /home/rpb/omicron/dataset-project-polymers/melt-v048-dox/a02-aside/source.gro
"""},
}
