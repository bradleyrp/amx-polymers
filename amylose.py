#!/usr/bin/env python

"""
Turn maltoheptaose into amylose.
"""

### STOLEN

### INTERLUDE to build a linked-list-style representation of the molecule to automatically find bonds

### automatically catalog all possible bonds (including angles and dihedrals)
class Molecule:
	def __init__(self):
		self.atoms = []
	def __repr__(self):
		return str(self.atoms)
	def get_atom(self,num):
		matches = [a for a in self.atoms if a.num==num]
		if len(matches)==0: raise Exception('cannot find atom %s'%num)
		elif len(matches)>1: raise Exception('dev. uniqueness violated!')
		else: return matches[0]
	def add_link(self,bond):
		i,j = [self.get_atom(b) for b in bond]
		if j not in i.neighbors: i.neighbors.append(j)
		if i not in j.neighbors: j.neighbors.append(i)
class Atom:
	def __init__(self,num):
		self.num = num
		self.neighbors = []
	def __repr__(self):
		return '%s: %s'%(self.num,[i.num for i in self.neighbors])

### END STOLEN

### BEGIN BACKWASH

import numpy as np
boxstuff = lambda pts,vec : pts-(pts>vec)*vec+(pts<np.array([0.,0.,0.]))*vec

def trim_waters(structure='solvate-dense',gro='solvate',gap=3,boxvecs=None,method='aamd',boxcut=True):
	"""
	ABSOLUTE FINAL VERSION OF THIS FUNCTION HOPEFULLY GOD WILLING
	trim_waters(structure='solvate-dense',gro='solvate',gap=3,boxvecs=None)
	Remove waters within a certain number of Angstroms of the protein.
	#### water and all (water and (same residue as water within 10 of not water))
	note that we vided the solvate.gro as a default so this can be used with any output gro file
	IS IT A PROBLEM THAT THIS DOESN'T TOUCH THE IONS??
	"""
	use_vmd = state.q('use_vmd',False)
	if (gap != 0.0 or boxcut) and use_vmd:
		if method == 'aamd': watersel = "water"
		elif method == 'cgmd': watersel = "resname %s"%state.q('sol')
		else: raise Exception("\n[ERROR] unclear method %s"%method)
		#---! gap should be conditional and excluded if zero
		vmdtrim = [
			'package require pbctools',
			'mol new %s.gro'%structure,
			'set sel [atomselect top \"(all not ('+\
			'%s and (same residue as %s and within '%(watersel,watersel)+str(gap)+\
			' of not %s)))'%watersel]
		#---box trimming is typical for e.g. atomstic protein simulations but discards anything outside
		if boxcut:
			vmdtrim += [' and '+\
			'same residue as (x>=0 and x<='+str(10*boxvecs[0])+\
			' and y>=0 and y<= '+str(10*boxvecs[1])+\
			' and z>=0 and z<= '+str(10*boxvecs[2])+')']
		vmdtrim += ['"]','$sel writepdb %s-vmd.pdb'%gro,'exit',]			
		with open(state.here+'script-vmd-trim.tcl','w') as fp:
			for line in vmdtrim: fp.write(line+'\n')
		vmdlog = open(state.here+'log-script-vmd-trim','w')
		#vmd_path = state.gmxpaths['vmd']
		vmd_path = 'vmd'

		#---previously used os.environ['VMDNOCUDA'] = "1" but this was causing segfaults on green
		p = subprocess.Popen('VMDNOCUDA=1 '+vmd_path+' -dispdev text -e script-vmd-trim.tcl',
			stdout=vmdlog,stderr=vmdlog,cwd=state.here,shell=True,executable='/bin/bash')
		p.communicate()
		#---!
		#with open(wordspace['bash_log'],'a') as fp:
		#	fp.write(gmxpaths['vmd']+' -dispdev text -e script-vmd-trim.tcl &> log-script-vmd-trim\n')
		gmx_run(state.gmxpaths['editconf']+' -f %s-vmd.pdb -o %s.gro -resnr 1'%(gro,gro),
			log='editconf-convert-vmd')
	#---scipy is more reliable than VMD
	elif gap != 0.0 or boxcut:
		import scipy
		import scipy.spatial
		import numpy as np
		#---if "sol" is not in the state we assume this is atomistic and use the standard "SOL"
		watersel = state.q('sol','SOL')
		incoming = read_gro(structure+'.gro')
		is_water = np.array(incoming['residue_names'])==watersel
		is_not_water = np.array(incoming['residue_names'])!=watersel
		water_inds = np.where(is_water)[0]
		not_water_inds = np.where(np.array(incoming['residue_names'])!=watersel)[0]
		points = np.array(incoming['points'])
		residue_indices = np.array(incoming['residue_indices'])
		if gap>0:
			#---previous method used clumsy/slow cdist (removed)
			#---use scipy KDTree to find atom names inside the gap
			#---note that order matters: we wish to find waters too close to not_waters
			#! note that cKDTree is OOM faster than KDTree
			# added boxsize which requires a recent-vintage scipy but includes PBCs in case out of box
			#! somewhat clumsy handling of box size
			#! switch to GMXStructure?
			boxsize = [float(i) for i in incoming['lines'][-1].split()]
			# put everything back in the box
			close_dists,neighbors = scipy.spatial.cKDTree(boxstuff(pts=points[not_water_inds],vec=boxsize),boxsize=boxsize).query(points[water_inds],distance_upper_bound=gap/10.0)
			#---use the distances to find the residue indices for waters that are too close 
			excludes = np.array(incoming['residue_indices'])[is_water][np.where(close_dists<=gap/10.0)[0]]
			#---get residues that are water and in the exclude list
			#---note that the following step might be slow
			exclude_res = [ii for ii,i in enumerate(incoming['residue_indices']) if i in excludes and is_water[ii]]
			#---copy the array that marks the waters
			surviving_water = np.array(is_water)
			#---remove waters that are on the exclude list
			surviving_water[exclude_res] = False
		else: 
			excludes = np.array([])
			surviving_water = np.ones(len(residue_indices)).astype(bool)
		#---we must remove waters that lie outside the box if there is a boxcut
		insiders = np.ones(len(points)).astype(bool)
		if boxcut:
			#---remove waters that lie outside the box
			#---get points that are outside of the box
			outsiders = np.any([np.any((points[:,ii]<0,points[:,ii]>i),axis=0) 
				for ii,i in enumerate(boxvecs)],axis=0)
			#---get residue numbers for the outsiders
			outsiders_res = np.array(incoming['residue_indices'])[np.where(outsiders)[0]]
			#---note that this is consonant with the close-water exclude step above (and also may be slow)
			exclude_outsider_res = [ii for ii,i in 
				enumerate(incoming['residue_indices']) if i in outsiders_res]
			insiders[exclude_outsider_res] = False
		surviving_indices = np.any((is_not_water,np.all((surviving_water,insiders),axis=0)),axis=0)
		lines = incoming['lines']
		lines = lines[:2]+list(np.array(incoming['lines'][2:-1])[surviving_indices])+lines[-1:]
		xyzs = list(points[surviving_indices])
		write_gro(lines=lines,xyzs=xyzs,output_file=state.here+'%s.gro'%gro)
	else: raise Exception('you need to either trim the box or remove waters in a gap')

#! note that this was fixed to use ckDtree and PBCs but the box is still too big
def solvate(structure,gro,edges=None,center=False):
	"""
	Standard solvate procedure for atomistic protein in water.
	Often requires restuff above.
	"""
	#---purge the wordspace of solvent and anions in case we are resuming
	for key in [state.q('anion'),state.q('cation'),'SOL']:
		if key in list(zip(*state.composition))[0]:
			del state.composition[list(zip(*state.composition))[0].index(key)]
	#---make an oversized water box
	boxdims_old,boxdims = get_box_vectors(structure)
	newdims = boxdims_old
	#---impose the total box thickness here!
	if state.thickness and edges: raise Exception('cannot set both thickness and edges')
	if state.q('thickness'): 
		boxdims[2] = state.q('thickness')
		boxdims_old[2] = state.q('thickness')
	if edges:
		#---custom water buffers in each direction require that we make a new box
		water_edges = [float(j) for j in state.water_edges]
		if not len(water_edges)==3: raise Exception('water_edges must be a triplet')
		boxdims_old = [boxdims_old[jj]+2*water_edges[jj] for jj in range(3)]
		structure_new = structure+'-centered'
		gmx('editconf',structure=structure,gro=structure_new,
			box=' '.join([str(j) for j in boxdims_old]),c=True,log='editconf-recenter')
		structure = structure_new
	#---if no solvent argument we go and fetch it
	solvent = state.q('solvent','spc216')
	if solvent=='spc216' and not os.path.isfile(state.here+'spc216.gro'):
		share_dn = gmx_get_share()
		shutil.copyfile(os.path.join(share_dn,'spc216.gro'),state.here+'spc216.gro')
	#---! solvent must be centered. for some reason spc216 is not in the box entirely.
	if solvent=='spc216':
		_,boxdims_spc216 = get_box_vectors('spc216')
		gmx('editconf',structure='spc216',gro='spc216-center',
			center=' '.join([str(i/2.) for i in boxdims_spc216]),log='spc216-recenter')
		solvent = 'spc216-center'
	#---get the basedim for the incoming water box
	basedim,_ = get_box_vectors(solvent)
	#---use the preexisting box vectors
	#---! fixed this from newdims to boxdims_old since the solvate function works in-place
	nbox = ' '.join([str(int(i/basedim[ii]+1)) for ii,i in enumerate(boxdims_old)])
	gmx('genconf',structure=solvent,gro='solvate-empty-uncentered-untrimmed',nbox=nbox,log='genconf')
	gro_combinator(structure+'.gro','solvate-empty-uncentered-untrimmed.gro',
		box=boxdims_old,cwd=state.here,gro='solvate-dense')
	atom_resolution = atomistic_or_coarse()
	trim_waters(structure='solvate-dense',gro=gro,gap=state.q('water_buffer',3),
		boxvecs=boxdims_old,method=atom_resolution,boxcut=True)
	#---! ugly
	sol = state.q('sol','SOL')
	nwaters = count_molecules(gro,sol)/({'aamd':3.0,'cgmd':1.0}[atom_resolution])
	if round(nwaters)!=nwaters: raise Exception('[ERROR] fractional water molecules')
	else: nwaters = int(nwaters)
	component(sol,count=nwaters)
	state.bilayer_dimensions_solvate = boxdims_old
	state.water_without_ions = nwaters

#! version for melts
def solvate_weird(structure,gro):
	"""
	"""
	# purge the wordspace of solvent and anions in case we are resuming
	for key in [state.q('anion'),state.q('cation'),'SOL']:
		if key in list(zip(*state.composition))[0]:
			del state.composition[list(zip(*state.composition))[0].index(key)]
	# make an oversized water box
	boxdims_old,boxdims = get_box_vectors(structure)
	newdims = boxdims_old
	# if no solvent argument we go and fetch it
	solvent = state.q('solvent','spc216')
	if solvent=='spc216' and not os.path.isfile(state.here+'spc216.gro'):
		share_dn = gmx_get_share()
		shutil.copyfile(os.path.join(share_dn,'spc216.gro'),state.here+'spc216.gro')
	#! solvent must be centered. for some reason spc216 is not in the box entirely.
	if solvent=='spc216':
		_,boxdims_spc216 = get_box_vectors('spc216')
		gmx('editconf',structure='spc216',gro='spc216-center',
			center=' '.join([str(i/2.) for i in boxdims_spc216]),log='spc216-recenter')
		solvent = 'spc216-center'
	# get the basedim for the incoming water box
	basedim,_ = get_box_vectors(solvent)
	# use the preexisting box vectors
	#! fixed this from newdims to boxdims_old since the solvate function works in-place
	nbox = ' '.join([str(int(i/basedim[ii]+1)) for ii,i in enumerate(boxdims_old)])
	sol = state.q('sol','SOL') #! ugly
	#! using solvate instead of genconf because it has maxsol but it's actually really basic. it just stacks
	#! ... up the boxes along edge. right now we have no choice because the molecules are long
	gmx('solvate',structure=structure,solvent=solvent,gro=gro,
		maxsol=int(component('AmylBig')*float(state.water_ratio)*state.n_p*state.beads_per_monomer),
		log='solvate')
	atom_resolution = atomistic_or_coarse()
	nwaters = count_molecules(gro,sol)/({'aamd':3.0,'cgmd':1.0}[atom_resolution])
	if round(nwaters)!=nwaters: raise Exception('[ERROR] fractional water molecules')
	else: nwaters = int(nwaters)
	component(sol,count=nwaters)
	state.bilayer_dimensions_solvate = boxdims_old
	state.water_without_ions = nwaters

### END BACKWASH

def make_amylose():
	"""
	Build amylose from our interpretation of maltoheptaose.
	See inputs/martini/martini-extensions/martini_v2.0_sugars.itp
	"""
	n_p = settings.n_p
	source = GMXTopology(settings.heptaose_source)
	from copy import deepcopy
	mol = deepcopy(source.molecules['Maltoheptaose'])
	if not (n_p%2)==1 or n_p<=7: raise Exception('n_p must be larger than 7 and odd: %d'%n_p)
	#! custom repetitor based on pattern in the maltoheptaose
	atoms = deepcopy(mol['atoms'][:3])+[i for j in 
		[deepcopy(mol['atoms'][3:9]) for k in range(((n_p-1)/2))] for i in j]
	#! hardcoding monomer beads here for a check in a few lines
	n_monomer_beads = 3
	# fix numbering
	atom_id,cgnr = 0,0
	for mono in range(n_p):
		for a in range(3):
			atoms[mono*3+a]['id'] = mono*3+a+1
			atoms[mono*3+a]['cgnr'] = mono+1
	# create abstract bonds
	bonds_abstract = []
	for mono in range(n_p):
		for i,j in [[1,2],[1,3],[1,4]]:
			this_bond = {'i':mono*3+i,'j':mono*3+j}
			if any([k>n_p*n_monomer_beads for k in this_bond.values()]): continue
			else: bonds_abstract.append({'i':mono*3+i,'j':mono*3+j})
	# MORE STOLEN
	molecule = Molecule()
	bond_list = [tuple([bond[k] for k in 'ij']) for bond in bonds_abstract]
	for atom_id in list(set([i for j in bond_list for i in j])):
		molecule.atoms.append(Atom(atom_id))
	for bond in bond_list: molecule.add_link(bond)
	def get_paths(atoms,path=None,d=3):
		if d>0:
			if not path: path,next_up = [],atoms
			else: next_up = path[-1].neighbors
			for neighbor in next_up:
				local_path = path[:]+[neighbor]
				for p in get_paths(atoms,path=local_path,d=d-1): yield p
		else: yield path
	# END MORE STOLEN
	do_constraints = False
	fkey = 'fc' if do_constraints else 'force'
	bond_lists = {
		'bonds':{
			('B2','B3'):{'funct':1,'length':0.222,fkey:30000},
			('B2','B1'):{'funct':1,'length':0.246,fkey:30000},
			('B2','B4'):{'funct':1,'length':0.561,fkey:30000},
			('B4','B6'):{'funct':1,'length':0.239,fkey:30000},
			('B4','B5'):{'funct':1,'length':0.281,fkey:30000},
			('B4','B2'):{'funct':1,'length':0.561,fkey:30000},
			},
		'angles':{
			('B3','B2','B4'):{'funct':2,'angle':150.,'force':50.},
			('B1','B2','B4'):{'funct':2,'angle':150.,'force':50.},
			('B2','B4','B6'):{'funct':2,'angle':150.,'force':50.},
			('B2','B4','B5'):{'funct':2,'angle':150.,'force':50.},
			},
		'dihedrals':{ # mult is 1
			('B6','B4','B2','B3'):{'funct':1,'angle':20.,'force':15.,'multiplicity':1},
			('B6','B4','B2','B1'):{'funct':1,'angle':55.,'force':15.,'multiplicity':1},
			('B5','B4','B2','B3'):{'funct':1,'angle':42.,'force':15.,'multiplicity':1},
			},
		}
	bond_types = ['bonds','angles','dihedrals']
	bonds = dict([(key,[]) for key in bond_types])
	for n_atoms,bond_type in zip([2,3,4],bond_types):
		bond_list = bond_lists[bond_type]
		for bond in get_paths(molecule.atoms,d=n_atoms):
			# get atom names for the bond
			this_inds = tuple([[int(a['id']) for a in atoms].index(i.num) for i in bond])
			this_bond_names = tuple([atoms[ind]['atom'] for ind in this_inds])
			# lookup the bond
			if this_bond_names in bond_list: this_bond = deepcopy(bond_list[this_bond_names])
			elif this_bond_names[::-1] in bond_list: this_bond = deepcopy(bond_list[this_bond_names[::-1]])
			else: 
				if bond_type=='bonds':
					print(bond)
				continue #! continue not pass! or everything is confused!
			these = [atoms[k] for k in this_inds]
			detail = dict([(k,these[kk]['id']) for kk,k in enumerate('ijkl'[:n_atoms])])
			this_bond.update(**detail)
			#if bond_type=='angle' and set(this_bond.keys())=={'i','length','fc','j'}
			bonds[bond_type].append(this_bond)

	#! hacking through some issues here!
	# remove redundancies
	#! note that this is corrects a design flaw from the way I have coded everything up above
	for bond_type in bond_types:
		bond_set_nonredundant,mirror = [],[]
		inds = [[i[k] for k in ('ijkl')[:dict(zip(bond_types,[2,3,4]))[bond_type]]] 
			for i in bonds[bond_type]]
		for ind,bond in zip(inds,bonds[bond_type]):
			keys = ('ijkl')[:dict(zip(bond_types,[2,3,4]))[bond_type]]
			if ind not in mirror and ind[::-1] not in mirror: 
				mirror.append(ind)
				bond_set_nonredundant.append(bond)
		bonds[bond_type] = bond_set_nonredundant

	mol.pop('bonds',None)
	mol['atoms'] = atoms
	mol['bonds'] = bonds['bonds']
	mol['angles'] = bonds['angles']
	mol['dihedrals'] = bonds['dihedrals']
	mol['moleculetype']['molname'] = settings.molecule_name
	new_top = GMXTopology()
	#import ipdb;ipdb.set_trace()
	new_top.add_molecule(**{settings.molecule_name:mol})
	new_top.write(state.here+'amylose31.itp')
	if 'itp' not in state: state.itp = []
	state.itp.append('amylose31.itp')
