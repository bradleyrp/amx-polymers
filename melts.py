#!/usr/bin/env python

"""
Utility functions for making polymer melts.
Oriented around the problem of getting accurate parameters for dextran, an oligoglucose 
composed of 1,6-glycosidic linkages not yet standard in MARTINI.
"""

import os,re,subprocess
import numpy as np
import MDAnalysis
from melts_simple import make_off_lattice_walk
from melts_tools import vecnorm,dotplace,rotation_matrix,align,apply_rotation,review3d
from melts_tools import planeproject,vecangle,measure_torsions

_not_reported = ['vecnorm','dotplace','rotation_matrix',
	'review3d','apply_rotation','align','vecangle','measure_torsions',
	'planeproject','vecangle','measure_torsions']

class MultiScaleModel:
	"""
	Manage the mappings between scales.
	"""
	def __init__(self,**kwargs):
		"""Unpack a mapping into a model."""
		mapping = kwargs.pop('mapping')
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		style = mapping.pop('style',None)
		# mappings which supply different "parts" (i.e. terminus or middle) of a polymer
		if style == 'coarse_to_atomistic': 
			self.style = style
			self.parts = mapping.pop('by_part',{})
		else: raise Exception('invalid style %s'%style)
		if mapping: raise Exception('unprocessed mapping directives %s'%mapping)
		# shared variables
		self.reps = {}

	def add_representation(self,name,**kwargs):
		"""Add a representation to the model."""
		if name in self.reps: raise Exception('representation %s already exists'%name)
		self.reps[name] = kwargs

	def build_fine_to_coarse(self,rep,target_name,**mods):
		"""
		Given a fine-grained representation and a mapping, prepare the coarse model and the indexer.
		"""
		source = self.reps[rep]
		if target_name in self.reps: raise Exception('')
		# molecule map can be overridden
		#! note that once we do two calls to build_fine_to_coarse for "bigger" we can make the 
		#! ... if statement below strict again and remove the base
		molecule_map = mods.get('molecule_map',source['molecule_map'])
		connector = lambda resnum,total: {0:None,1:'terminus_start',
			total:'terminus_end'}.get(resnum,'middle')
		#---the assignment of "parts" descriptors depends on the mapping type
		if molecule_map in ['one molecule many residues injective',
			'one molecule many residues bigger']:
			"""
			In this section we perform the residue-injective mapping.
			We build a list of beads with names and residues taken from the unique atomistic residues.
			Each coarse residue has a list of atomistic indices *from that residue index* that map to it.
			"""
			#---get unique resnums because this mapping is injective between residues
			source['resnums_u'] = np.unique(source['resnums'])
			#---use the molecule map style to assign parts to each atom
			source['parts'] = [connector(r,len(source['resnums_u'])) for r in source['resnums']]
			#---build the coarse model
			target = dict(resnums_u=source['resnums_u'],beads=[])
			#---loop over resnums/coarse-grained residues
			for resnum in target['resnums_u']:
				#---each bead gets a residue index
				bead = dict(resnum=resnum)
				#---determine the part
				parts = [source['parts'][i] for ii,i in enumerate(np.where(source['resnums']==resnum)[0])]
				parts = list(set(parts))
				if len(parts)!=1: 
					raise Exception('in the injective residue method we failed to assign a unique part')
				else: bead['part'] = parts[0]
				#---get the beads for the part
				self.parts[bead['part']]
				#---loop over index pairs
				for names_coarse,names_fine in self.parts[bead['part']]:
					#---currently this allows only coarse-to-fine
					if len(names_coarse)!=1: raise Exception('development')
					specify = dict(name=names_coarse[0])
					specify['indices'] = np.where(np.all((np.in1d(source['names'],names_fine),
						source['resnums']==resnum),axis=0))[0]
					target['beads'].append(dict(bead,**specify))
				#---having assembled a list of beads we repackage them to mimic the incoming data
				target['resnums'] = [i['resnum'] for i in target['beads']]
				target['names'] = [i['name'] for i in target['beads']]
			#---make the final mapping
			target['coarsen'] = [i['indices'] for i in target['beads']]
			target['length'] = len(target['coarsen'])
		# very similar to the above block with no mapping because there is none from a pentamer
		# ... to a much larger structure. bond distributions are generalized later.
		# we actually overwrite this after getting the coarsening above
		if molecule_map == 'one molecule many residues bigger':
			nres = mods['nres']
			# write the original resnames
			target['resnums_u_base'] = target['resnums_u']
			target['resnums_u'] = np.arange(nres)+1
			target['beads_base'] = target['beads']
			target['beads'] = []
			# loop over resnums/coarse-grained residues
			for resnum in target['resnums_u']:
				# each bead gets a residue index
				bead = dict(resnum=resnum,part=connector(resnum,nres))
				# get the beads for the part
				self.parts[bead['part']]
				# loop over index pairs
				for names_coarse,names_fine in self.parts[bead['part']]:
					# currently this allows only coarse-to-fine
					if len(names_coarse)!=1: raise Exception('development')
					specify = dict(name=names_coarse[0])
					specify['indices'] = np.where(np.all((np.in1d(source['names'],names_fine),
						source['resnums']==resnum),axis=0))[0]
					target['beads'].append(dict(bead,**specify))
		target['molecule_map'] = molecule_map
		self.reps[target_name] = target

	def reduce(self,source,target,positions):
		"""Convert a set of atomistic positions to coarse ones."""
		pos = [np.average(positions[i],axis=0,weights=self.reps[source]['masses'][i]) 
			for i in self.reps[target]['coarsen']]
		return pos

	def interpret_model(self,source,coords,**model_spec):
		"""
		Read an abstract topology definition.
		We will measure features of the model from the coords and use it to construct the topology.
		"""
		rep = self.reps[source]
		molecule_map = rep['molecule_map']
		remember_bonds = False
		# note that we automatically choose general bonds instead of explicit ones
		if molecule_map=='one molecule many residues injective':
			general_bonds = False
			remember_bonds = True
		elif molecule_map=='one molecule many residues bigger':
			general_bonds = True
			#! hard-coded reference for now
			rep_ref = self.reps['coarse']
		# do not alter general_bonds above, since this logic must be strict!
		else: raise Exception
		remember_bonds = model_spec.pop('remember_bonds',remember_bonds)
		# collect the intra-residue bonds
		bonds_intra,bond_parts = [],[]
		#! this is currently hard-coded. check the style later
		model_by_part = model_spec['model']['by_part']
		#! insert a check to make sure the model has all the parts
		#! save partnames to the datastructure?
		partnames = list(set([i['part'] for i in rep['beads']]))

		# make a blank GMXTopology (note that you must always include the name here)
		new_top = GMXTopology()

		# per-atom topology column requirements to be read from the model spec
		requirements = {
			'atoms':['charge','mass','type'],}
		# quicker lookups between residues
		def get_bead_in_residue(bead_name,resnum,beads_name='beads'):
			matches = [ii for ii,i in enumerate(rep[beads_name]) 
				if i['name']==bead_name and i['resnum']==resnum]
			if len(matches)>1: raise Exception
			elif len(matches)==0: return None
			else: return matches[0]
		# construct the atoms list which is required for indexing the bonds later on
		atom_id,atoms = 1,[]
		bonds,angles,dihedrals,bonds_abstract = [],[],[],[]
		# loop over residues
		for resnum in rep['resnums_u']:
			# get the beads in this residue
			beads_inds,beads = zip(*[(ii,i) for ii,i in enumerate(rep['beads']) if i['resnum']==resnum])
			#! first bead is the resname
			resname = beads[0]['part']
			# loop over beads
			for atom_residue_index,(index,bead) in enumerate(zip(beads_inds,beads)):
				# construct an atom entry
				# several columns are automatic or identical to the abstract model
				new_atom = dict(id=atom_id,resnr=resnum,cgnr=resnum)
				atom_id += 1
				#! is this elegant?
				new_atom.update(resname=model_spec['model']['by_part'][resname]['resname'])
				new_atom.update(atom=bead['name'])
				# collect atom-specific columns from the model spec
				for key in requirements['atoms']:
					#! ryan is worried this is still too complicated?
					new_atom[key] = model_spec['model']['by_part'][resname][
						'atoms'][atom_residue_index][bead['name']][key]
				atoms.append(new_atom)
			# get connectivity here
			for bead_1,bead_2 in model_spec['model']['by_part'][resname].get('bonds',[]):
				lmatch = get_bead_in_residue(bead_name=bead_1,resnum=resnum)
				rmatch = get_bead_in_residue(bead_name=bead_2,resnum=resnum)
				bead_1_i,bead_2_i = lmatch,rmatch
				new_bond = dict(i=bead_1_i+1,j=bead_2_i+1)			
				bonds_abstract.append(new_bond)
		# handle bonds between residues
		resnums_adjacent = [(rep['resnums_u'][i],rep['resnums_u'][i+1]) 
			for i in range(len(rep['resnums_u'])-1)]
		between_dict = model_spec['model'].get('between_parts',{})
		for left,lspec in between_dict.items():
			for right,rspec in lspec.items():
				# search for all bonds between adjacent residues
				for bead_1,bead_2 in rspec.get('bonds',[]):
					# check all adjacent residues
					for r1,r2 in resnums_adjacent:
						lmatch = get_bead_in_residue(bead_name=bead_1,resnum=r1)
						rmatch = get_bead_in_residue(bead_name=bead_2,resnum=r2)
						if lmatch and rmatch: 
							bead_1_i,bead_2_i = lmatch,rmatch
							new_bond = dict(i=bead_1_i+1,j=bead_2_i+1)
							bonds_abstract.append(new_bond)

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
				i,j = [molecule.get_atom(b) for b in bond]
				if j not in i.neighbors: i.neighbors.append(j)
				if i not in j.neighbors: j.neighbors.append(i)
		class Atom:
			def __init__(self,num):
				self.num = num
				self.neighbors = []
			def __repr__(self):
				return '%s: %s'%(self.num,[i.num for i in self.neighbors])
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

		# save bond distributions for inspection
		bonds_review = dict([(i,{}) for i in [2,3,4]])
		# loop over angles and dihedrals
		for nbeads in [2,3,4]:
			bead_candidates = []
			for a in get_paths(molecule.atoms,d=nbeads):
				this = [i.num for i in a]
				if (this not in bead_candidates and this[::-1] not in bead_candidates 
					and len(set(this))==nbeads):
					bead_candidates.append(this)
			valid_bonds = bead_candidates
			# loop over these bonds
			for beads in sorted(valid_bonds):
				if nbeads==2:
					ii,jj = beads
					inds = [m-1 for m in beads]
					new_bond = dict(i=ii,j=jj)
					if not general_bonds:
						#---get the positions for this bond pair
						positions = coords[:,np.array(inds)]
						#---compute observed distances
						distances = np.linalg.norm(positions.ptp(axis=1),axis=1)
						bonds_review[2][(ii,jj)] = distances
						new_bond['length'] = distances.mean()
					else:
						# previously we got the exact index on coords (reduced i.e. coarse-grained positions) 
						# ... from the atomistic simulation. now we read the coords by name to generalize it
						# get all rows in coods that match the name of bead_1 and bead_2
						bead_1,bead_2 = [rep['beads'][b]['name'] for b in inds]
						lmatch,rmatch = [np.array([i for i in [get_bead_in_residue(
							bead_name=bead,resnum=r,beads_name='beads_base') 
							for r in rep.get('resnums_u_base',rep['resnums_u'])] if i!=None]) 
							for bead in [bead_1,bead_2]]
						distances = (coords[:,lmatch]-coords[:,rmatch])
						#! note that distances is nframes by the number of redudant beads of that type by XYZ
						distances_norm = np.array([np.linalg.norm(bd.T,axis=1).mean() 
							for bd in distances.T])
						new_bond['length'] = distances_norm.mean()
						bonds_review[2][(ii,jj)] = distances_norm
					new_bond['funct'] = 1
					new_bond['force'] = -1
					bonds.append(new_bond)
				elif nbeads==3:
					#! beware python 2 will reset i in a comprehension hence this is a poor index choice
					ii,jj,kk = beads
					#! back to zero indexes
					inds = [m-1 for m in beads]
					new_angle = dict(i=ii,j=jj,k=kk)
					new_angle['funct'] = 2
					#! started with 2000 which was unstable. then 500. then used common sense to check
					#! ... martini_v2.0_sugars.itp and found that most angles are type 2 with forces
					#! ... around 25, 50, 100. some go up to 800 for e.g. Laminaraheptabiose
					# replaced force here because it is tuned below
					new_angle['force'] = -1.
					# measure the angles 1-1 from the underlying coordinates
					if not general_bonds:
						vecs1,vecs2 = np.array((
							coords[:,inds[0]]-coords[:,inds[1]],
							coords[:,inds[2]]-coords[:,inds[0]]))
						angles_measured = np.array([vecangle(i,j) for i,j in zip(vecs1,vecs2)])
						bonds_review[3][(ii,jj,kk)] = angles_measured
						new_angle['angle'] = angles_measured.mean()
					# the coarse model may have more coordinates than the (reduced) source
					else: 
						# to get the generic angles we need to check the inds list against the valid bond list
						# ... because just checking by name means we end up looking for all of the SB2, then
						# ... all of the SB3, then all of the MB1 and we get 1,1,3 candidates. instead of this
						# ... method we need to generalize the entire angle at once.
						bead_names = [rep['beads'][b]['name'] for b in inds]
						#! vital note: we prepared to use resnums to ensure that the bead names were also 
						#! ... split across residue numbers correctly. however, this was not needed since the 
						#! ... following  section only selects triplets from the reference structure and 
						#! ... these are always split across residues correctly. basically the for/if loop 
						#! ... below that picks the right angles always ensures we get valid angles as long 
						#! ... as both models have them. this might not always be the case in branched 
						#! ... topologies
						resnums = [rep['beads'][b]['resnum'] for b in inds]
						# by the time you arrive here we have a molecule map "bigger" which requires
						# ... general bonds and we need to match the coarse_big rep against the coarse rep
						valid_angles_ref = [tuple([a[k] for k in 'ijk']) 
							for a in rep_ref['bond_info']['angles']]
						# search valid angles from the reference for those that match our general request
						angle_matches = []
						for angle in valid_angles_ref:
							beads_ref = [rep_ref['beads'][i-1] for i in angle]
							if ([b['name'] for b in beads_ref]==bead_names or 
								[b['name'] for b in beads_ref][::-1]==bead_names):
								angle_matches.append([i-1 for i in angle])
						matches = [np.array([i for i in [get_bead_in_residue(
							bead_name=bead,resnum=r,beads_name='beads_base') 
							for r in rep.get('resnums_u_base',rep['resnums_u'])] if i!=None]) 
							for bead in bead_names]
						# observations are instance of angle type by frames by bead by XYZ
						# ... note that for middle beads we have more than one instance of the angle type
						obs = np.array([coords[:,i] for i in np.array(angle_matches)])
						vecsAB = obs[:,:,2]-obs[:,:,1]
						vecsCB = obs[:,:,0]-obs[:,:,1]
						angles_measured = np.array([np.array([vecangle(i,j) 
							for i,j in zip(vecsAB[iii],vecsCB[iii])]) for iii in range(len(vecsAB))])
						bonds_review[3][(ii,jj,kk)] = angles_measured
						# average over instances and frames
						new_angle['angle'] = angles_measured.mean()
					angles.append(new_angle)
				elif nbeads==4:
					ii,jj,kk,ll = beads
					if not general_bonds:
						#! back to zero indexes
						inds = [m-1 for m in beads]
						new_dihedral = dict(i=ii,j=jj,k=kk,l=ll)
						new_dihedral['funct'] = 1
						new_dihedral['multiplicity'] = 1 # MAY NOT REALLY BE ONE!
						#! checked martini_v2.0_sugars.itp and saw fc was 5, 15, etc
						# force tuning how happens below
						new_dihedral['force'] = -1.
						dihedrals_measured = measure_torsions(coords,[inds])
						bonds_review[4][(ii,jj,kk,ll)] = dihedrals_measured
						new_dihedral['angle'] = dihedrals_measured.mean()
						dihedrals.append(new_dihedral)
					else:
						# mimic the angle section
						inds = [m-1 for m in beads]
						new_dihedral = dict(i=ii,j=jj,k=kk,l=ll)
						bead_names = [rep['beads'][b]['name'] for b in inds]
						valid_dihedrals_ref = [tuple([a[k] for k in 'ijkl']) 
							for a in rep_ref['bond_info']['dihedrals']]
						# search valid angles from the reference for those that match our general request
						dihedral_matches = []
						for dihedral in valid_dihedrals_ref:
							beads_ref = [rep_ref['beads'][i-1] for i in dihedral]
							if ([b['name'] for b in beads_ref]==bead_names or 
								[b['name'] for b in beads_ref][::-1]==bead_names):
								dihedral_matches.append([i-1 for i in dihedral])
						dihedrals_measured = measure_torsions(coords,dihedral_matches)
						new_dihedral['funct'] = 1
						new_dihedral['multiplicity'] = 1 # MAY NOT REALLY BE ONE!
						new_dihedral['force'] = 5.
						bonds_review[4][(ii,jj,kk,ll)] = dihedrals_measured
						new_dihedral['angle'] = dihedrals_measured.mean()
						dihedrals.append(new_dihedral)
				else: raise Exception('development note. cannot generate abond with %d beads'%nbeads)

		# save the bond distributions
		import scipy
		import scipy.stats
		bond_stats = dict([(i,{}) for i in [2,3,4]])
		for anum,nbeads in enumerate([2,3,4]):
			vals_cat = np.concatenate(bonds_review[nbeads].values())
			bins = np.linspace(vals_cat.min(),vals_cat.max(),100)
			for bnum,bond in enumerate(bonds_review[nbeads].keys()):
				obs = bonds_review[nbeads][bond]
				counts,edges = np.histogram(obs,bins=bins,normed=True)
				mids = (edges[1:]+edges[:-1])/2.
				mu,std = scipy.stats.norm.fit(obs)
				bond_stats[nbeads][bond] = dict(mu=mu,std=std)
		if settings.get('do_review_plots',False):
			import matplotlib as mpl
			import matplotlib.pyplot as plt
			from plotter_omni_panels import square_tiles
			n_bonds = sum([len(bonds_review[i].keys()) for i in [2,3,4]])
			# report the bond distributions
			axes,fig = square_tiles(n_bonds,figsize=(15.),favor_rows=True,wspace=0.4,hspace=0.4)
			ax_count = 0
			for anum,nbeads in enumerate([2,3,4]):
				vals_cat = np.concatenate(bonds_review[nbeads].values())
				bins = np.linspace(vals_cat.min(),vals_cat.max(),100)
				for bnum,bond in enumerate(bonds_review[nbeads].keys()):
					obs = bonds_review[nbeads][bond]
					counts,edges = np.histogram(obs,bins=bins,normed=True)
					mids = (edges[1:]+edges[:-1])/2.
					mu,std = scipy.stats.norm.fit(obs)
					ax = axes[ax_count]
					ax.plot(mids,counts,c='k')
					ax.plot(mids,scipy.stats.norm.pdf(mids,mu,std),c='r',zorder=2)
					ax.set_title({2:'bond',3:'angle',4:'dihedral'}[nbeads]+' '+'-'.join([str(i) 
						for i in bond]))
					ax_count += 1
			# plot histograms of the widths
			fig = plt.figure()
			axes = [plt.subplot(k) for k in [211,212]]
			for anum,nbeads in enumerate([3,4]):
				vals = [v['std'] for k,v in bond_stats[nbeads].items()]
				counts,edges = np.histogram(vals,bins=len(vals))
				ax = axes[anum]
				ax.bar((edges[1:]+edges[:-1])/2.,counts)
				ax.set_title('%s std'%{2:'bond',3:'angle',4:'dihedral'}[nbeads])
				ax.set_xlabel('std')
				ax.set_ylabel('counts')
			plt.savefig('bonds_review_widths%s.png'%('_general' if general_bonds else ''))



		from makeface import import_remote
		incoming_filter_functions = import_remote(settings.bond_tuners_code)
		filter_bonds = incoming_filter_functions['filter_bonds']
		filter_angles = incoming_filter_functions['filter_angles']
		filter_dihedrals = incoming_filter_functions['filter_dihedrals']
		scale_bonds = incoming_filter_functions.get('scale_bonds',lambda x:x)

		# only filter the bonds if we are not remembering them for a second pass
		if not remember_bonds:
			# filter over bond types
			for nbeads in [2,3,4]:
				reworked = []
				subjects = {2:bonds,3:angles,4:dihedrals}[nbeads]
				filter_func = {2:filter_bonds,3:filter_angles,4:filter_dihedrals}[nbeads]
				letters = {2:'ij',3:'ijk',4:'ijkl'}[nbeads]
				for subject in subjects:
					names = [[a['atom'] for a in atoms if a['id']==subject[l]][0] for l in letters]
					resids = [[a['resnr'] for a in atoms if a['id']==subject[l]][0] for l in letters]
					strength = filter_func(bond_stats[nbeads][tuple([subject[k] 
						for k in letters])]['std'],names=names,resids=resids)
					if strength>0:
						reworked.append(subject)
						reworked[-1]['force'] = strength
						if nbeads==2:
							reworked[-1]['length'] = scale_bonds(reworked[-1]['length'],
								names=names,resids=resids,style='standard')
				print('[NOTE] %d/%d %s are remaining after the filter'%(
					len(reworked),len(subjects),{2:'bonds',3:'angles',4:'dihedrals'}[nbeads]))
				if nbeads==2: bonds = reworked
				elif nbeads==3: angles = reworked
				elif nbeads==4: dihedrals = reworked
				else: raise Exception

		#! eventually code an option for restraints here
		bonds_constraints = dict(bonds=bonds)
		if remember_bonds:
			self.reps[source]['bond_info'] = dict(angles=angles,dihedrals=dihedrals,**bonds_constraints)
		else:
			# populate the topology
			#! note that you have to pass in moleculetype or the ITP is incomplete
			new_top.molecules['DEX'] = dict(moleculetype=dict(molname='DEX',nrexcl=3),
				atoms=atoms,angles=angles,dihedrals=dihedrals,**bonds_constraints)
			if settings.get('no_terminals',False):
				mol = new_top.molecules['DEX']
				discards = range(1,3+1)+range(len(mol['atoms']),len(mol['atoms'])-3,-1)[::-1]
				new_atoms = []
				for ii in range(3,len(mol['atoms'])-3):
					new_atoms.append(mol['atoms'][ii])
					new_atoms[-1]['id'] -= 3
					new_atoms[-1]['resnr'] -= 1
				mol['atoms'] = new_atoms
				for btype in ['bonds','angles','dihedrals']:
					new_set = []
					letters = {'bonds':'ij','angles':'ijk','dihedrals':'ijkl'}
					if btype not in mol: continue
					for item in mol[btype]:
						if any([item[l] in discards for l in letters[btype]]): continue
						for l in letters[btype]: item[l] -= 3
						new_set.append(item)
					mol[btype] = new_set
			#! constraints
			new_top.write(state.here+'dextran.itp')
			#! writing multiple topologies
			mol['constraints'] = mol.pop('bonds')
			for c in mol['constraints']: c.pop('force')
			new_top.write(state.here+'dextran_constraints.itp')

### CONSTRUCTION FUNCTIONS (simple functions which build a starting structure)

def make_polymer(name='melt',**kwargs):
	"""
	Make a melt by placing monomers according to an off-lattice walk.
	This is the first section of make_cg_gel; the second section built the topology.
	"""
	polymer_molname = settings.molecule_name
	monomer_resname = settings.residue_name
	no_terminals = settings.get('no_terminals',False)
	# destination
	cwd = state.here
	n_p = state.melt_settings['n_p']
	angle,torsion = state.melt_settings['angle'],state.melt_settings['torsion']
	a0 = state.melt_settings['a0']
	# prepare an abstract 3D walk
	walk_abstract_pts = make_off_lattice_walk(n_p-1+(2 if no_terminals else 0),a0,angle,torsion)
	points_with_sides = []
	for pt in walk_abstract_pts:
		points_with_sides.append(pt)
		for i in range(2):
			#! randomly place sidechain points 1 a0 distance away. map 0 to 1 to -1 to 1
			points_with_sides.append(pt+(1-2*np.random.rand(3))/10.)
	points_with_sides = np.array(points_with_sides)
	residue_names = np.array([monomer_resname for p in points_with_sides])
	residue_indices = np.array([i/3+1 for i in range((n_p)*3)])
	#! hardcoding the naming scheme here but it would be useful to get this from the YAML file!
	#! working on no terminals
	if settings.get('no_terminals',False): 
		atom_names = np.array(['%sB%d'%(l,i) for l in 'S'+'M'*(n_p)+'E' for i in range(1,3+1)])
	else: atom_names = np.array(['%sB%d'%(l,i) for l in 'M'+'M'*(n_p)+'M' for i in range(1,3+1)])
	if not no_terminals:
		polymer = GMXStructure(pts=points_with_sides,residue_indices=residue_indices,
			residue_names=residue_names,atom_names=atom_names,box=[10.,10.,10.])
	else: 
		sl = slice(3,-3)
		polymer = GMXStructure(pts=points_with_sides[sl],residue_indices=residue_indices[sl],
			residue_names=residue_names[sl],atom_names=atom_names[sl],box=[10.,10.,10.])
	polymer.write(state.here+'%s-built.gro'%name)
	gmx('editconf',structure='%s-built'%name,gro=name,c=True,log='editconf-center-polymer')
	# add to the composition
	state.itp = ['dextran.itp']
	component(polymer_molname,count=1)

def make_single_polymer():
	"""
	Make very crude estimate for coordinates for a single polymer for minimization.
	To use this you must supply the correct ITP file and add it to the composition from your script.
	This method was originally designed for using stock MARTINI maltoheptaose.
	"""
	n_p = settings.n_p
	a0 = settings.a0
	bpm = settings.beads_per_monomer
	backbone_pos = settings.backbone_position
	residue_name = settings.residue_name
	#! note that we could send lambda functions directly in through YAML but they cannot write to state
	atom_namer = eval(settings.atom_namer)
	# backbone beads in a straight line along the x-direction
	backbone_coords = np.transpose([np.arange(n_p)*a0,np.zeros(n_p),np.zeros(n_p)])
	sidechains = np.array([[(1-2*np.random.rand(3))/10.+b for i in range(bpm-1)] for b in backbone_coords])
	# interleave the backbone and sidechain by backbone position
	bead_order = np.zeros(bpm)
	bead_order[backbone_pos-1] = 1
	coords = np.zeros((bpm*n_p,3))
	side_place = np.where(bead_order==0)[0]
	back_place = np.array([backbone_pos-1])
	for mnum in range(n_p):
		mono_coords_this = np.zeros((bpm,3))
		mono_coords_this[side_place] = sidechains[mnum]
		mono_coords_this[back_place] = backbone_coords[mnum]
		coords[mnum*bpm:(mnum+1)*bpm] = mono_coords_this
	atom_names = [atom_namer(i) for i in range(n_p*bpm)]
	residue_indices = np.ones(n_p*bpm).astype(int)
	residue_names = [residue_name for i in range(n_p*bpm)]
	box_extra = 2. # make the box bigger
	box_vectors = np.array([n_p*a0*box_extra for i in range(3)])
	coords -= coords.mean(axis=0)
	coords += box_vectors/2.
	# protect against tiny boxes
	box_vectors = [max([10.,i]) for i in box_vectors]
	polymer = GMXStructure(pts=coords,residue_indices=residue_indices,
		residue_names=residue_names,atom_names=atom_names,box=box_vectors)
	polymer.write(state.here+'vacuum_crude.gro')

def make_crude_coarse_polymer(name='vacuum',diagonals=False,review=False,**kwargs):
	"""
	Settings in the experiment bring you here, or above to `make_polymer`.
	The `make_polymer` code above will make a single large molecule. 
	This method will make a melt from a square lattice. 
	See make_single_polymer for an alternative.
	UNDER DEVELOPMENT.
	"""
	monomer_resname = settings.residue_name
	polymer_molname = settings.molecule_name
	### MAKE A SQUARE LATTICE (recapitulated from make_gel_on_lattice)
	cwd = state.here
	# make a universe
	a0 = state.lattice_melt_settings['a0']
	sizer = state.lattice_melt_settings['sizer']
	n_p = state.lattice_melt_settings['n_p']
	volume_limit = state.lattice_melt_settings['volume_limit']
	#! assume square
	xd,yd,zd = a0,a0,a0
	xld,yld,zld = sizer,sizer,sizer
	data_discrete = np.array([[i,j,k] 
		for k in np.arange(0,zld,1) 
		for j in np.arange(0,yld,1) 
		for i in np.arange(0,xld,1)]).astype(int)
	# generate neighbor rules (with periodic boundary conditions)
	if not diagonals:
		neighbor_list = [(0+i,0+j,0+k) 
			for i in range(-1,2) for j in range(-1,2) 
			for k in range(-1,2) if not all([x==0 for x in [i,j,k]]) 
			and sum([x!=0 for x in [i,j,k]])==1]
	else: raise Exception
	neighbors = lambda x,y,z: [((x+i)%xld,(y+j)%yld,(z+k)%zld) for i,j,k in neighbor_list]
	# CONSTRUCT A non-intersecting walk in 3-space
	# track traffic (no intersections) and random walks through the space
	walkers,walks = [],np.zeros((xld,yld,zld)).astype(int)
	walk_id = 1
	while float(np.sum(walks==0))/np.product(walks.shape)>volume_limit:
		filled_up = np.product(walks.shape)-float(np.sum(walks==0))
		status('[CONSTRUCT] making polymers',i=filled_up,looplen=np.product(walks.shape)*(1.0-volume_limit))
		prop_filled = float(np.sum(walks>0))/np.product(walks.shape)
		if prop_filled>=volume_limit: 
			status('breaking construction loop because the space-filled is %.3f'%prop_filled)
			break
		# random starting point
		pos = tuple(np.transpose(np.where(walks==0))[np.random.randint(np.sum(walks==0))])
		# initiate
		size = 0
		walks[pos] = walk_id
		# loop
		this_walk = [pos]
		failure = False
		while size<(n_p-1):
			ways = list(filter(lambda n : walks[n]==0,neighbors(*pos)))
			if not ways:
				failure = True
				break
			pos = ways[np.random.randint(len(ways))]
			this_walk.append(pos)
			walks[pos] = walk_id
			size += 1
		if not failure:
			walkers.append(this_walk)
			walk_id += 1
	# use a half-lattice offset to place the dots in the center of the box if we are doing a 3D review
	offset = a0/2.0
	# now that we have points, we recapitulate a sequence from make_cg_gel
	coords,atom_names,residue_indices = [],[],[]
	# loop over polymers
	for wnum,walker in enumerate(walkers):
		for pnum,pt in enumerate(walker):
			coords.append(pt)	
			#! coarse-grained mapping: three points per monomer
			for i in range(2):
				#! randomly place sidechain points 1 a0 distance away. map 0 to 1 to -1 to 1
				coords.append(pt+(1-2*np.random.rand(3))/10.)
				#! let's eliminate terminals for now!
				if not settings.get('no_terminals',False):
					if pnum==0: atom_name_letter = 'S'
					elif pnum<n_p-1: atom_name_letter = 'M'
					elif pnum==n_p-1: atom_name_letter = 'E'
					else: raise Exception
				else: atom_name_letter = 'M'
			atom_names.extend(['%sB%d'%(atom_name_letter,i) for i in range(1,3+1)])
			residue_indices.extend([wnum+1 for i in range(3)])
	coords = np.array(coords)*a0
	#! previously hacked the atom names for maltoheptaose but now we build with `gmx insert-molecules`
	atom_names = np.array(atom_names)
	residue_indices = np.array(residue_indices)
	residue_names = np.array([monomer_resname for p in coords])
	#! note that the box is dilated a bit to avoid a jammed system
	box_extra = 1.0
	polymer = GMXStructure(pts=coords,residue_indices=residue_indices,
		residue_names=residue_names,atom_names=atom_names,box=[sizer*a0*box_extra for i in range(3)])
	polymer.write(state.here+'%s.gro'%name)
	component(polymer_molname,count=len(walkers))
	#! also add the itp here
	#! moved to the script: state.itp = ['dextran.itp']

def hydration_adjust(structure,gro):
	"""
	Reduce the hydration level to a fraction of CG polymer beads to water beads.
	"""
	water_ratio = state.q('water_ratio',None)
	if water_ratio:
		sol_resname = settings.get('sol','SOL')
		polymer_molname = settings.molecule_name
		n_water_beads = dict(state.composition)[sol_resname]
		n_polymer_beads = float(state.n_p)*dict(state.composition)[polymer_molname]
		struct = GMXStructure(state.here+'%s.gro'%structure)
		remove_waters = n_water_beads - water_ratio*n_polymer_beads
		if remove_waters<0: raise Exception(('there are %s polymer beads and %s water beads so we cannot '
			'hydrate to a level of %s. try making the box bigger')%(
			n_polymer_beads,n_water_beads,water_ratio))
		water_inds = np.where(struct.residue_names==sol_resname)[0]
		n_keep_waters = int(water_ratio*n_polymer_beads)
		reindex = range(len(water_inds))
		np.random.shuffle(reindex)
		keeper_waters = np.sort(water_inds[reindex[:n_keep_waters]])
		keepers = np.sort(np.concatenate([np.where(struct.residue_names!=sol_resname)[0],keeper_waters]))
		struct.atom_names = struct.atom_names[keepers]
		struct.residue_indices = struct.residue_indices[keepers]
		struct.residue_names = struct.residue_names[keepers]
		struct.points = struct.points[keepers]
		struct.write(state.here+'%s.gro'%gro)
		component(sol_resname,count=n_keep_waters)
	else: copy_file('%s.gro'%structure,'%s.gro'%gro)

def write_structure_by_chain(structure='system',gro='system_chains'):
	"""
	Ensure that each chain has its own residue number. Useful for visualization.
	"""
	struct = GMXStructure(state.here+'%s.gro'%structure)
	polymer_molname = settings.molecule_name
	residue_name = settings.residue_name
	#! previously used DMR in the line below and DEX below that, in the component on the if not
	inds = np.where(struct.residue_names==residue_name)[0]
	#! assume that we only have a uniform melt of n_p, three beads per, minus two from the end points
	n_p = state.n_p
	bpm = settings.beads_per_monomer
	#! removed a -2 correction to the n_p for use in maltoheptaose which was necessary for removing terminals
	if not float(len(inds))/float(bpm)/(n_p)==component(polymer_molname): 
		raise Exception('failed to write structure by chains')
	resnums = (np.floor(np.arange(len(inds))/((n_p)*bpm))+1).astype(int)
	struct.residue_indices[inds] = resnums
	struct.write(state.here+'%s.gro'%gro)

### INTERFACE FUNCTIONS, called from the simulation script, which perform the mappings

def forward_mapper(write_coarse_coordinates=False,inspect_distributions=False):
	"""
	Compute statistics from an atomistic polymer on a coarse-grained mapping.
	Used to perform a "forward mapping" from atomistic to coarse-grained.
	"""
	# get the reference structure, presumably from another run
	ref_specs = settings.atomistic_reference
	ref_gro = os.path.join(ref_specs['path'],ref_specs['gro'])
	ref_xtc = os.path.join(ref_specs['path'],ref_specs['xtc'])
	uni = MDAnalysis.Universe(ref_gro,ref_xtc)
	sel = uni.select_atoms(ref_specs['selection'])
	
	# read the mapping
	import yaml
	with open(settings.mapping_spec) as fp: mapping = yaml.load(fp.read())
	model = MultiScaleModel(**mapping)

	#! move this into the class?
	lenscale = 10.0
	mass_table = {'H':1.008,'C':12.011,'O':15.999,'N':14.007,'P':30.974}
	# the atomistic reference tells us how to get molecules from that simulation
	molecule_map = ref_specs.get('molecule_map',None)
	# read the statistics for a single molecule
	# prepare a crude structure for the atomistic points
	#! this block might be included in the model class but for now it's only tangent
	fine = dict(molecule_map=molecule_map,selection=ref_specs['selection'])
	fine.update(resnums=sel.resnums,names=sel.names)
	fine.update(masses=np.array([mass_table[i[0]] for i in fine['names']]))
	model.add_representation('fine',**fine)
	nres = settings.melt_settings['n_p']+(2 if settings.get('no_terminals',False) else 0)
	# after preparing fine we make the coarse model
	if molecule_map == 'one molecule many residues injective':
		model.build_fine_to_coarse('fine','coarse')
	elif molecule_map == 'one molecule many residues bigger':
		model.build_fine_to_coarse('fine','coarse',molecule_map='one molecule many residues injective')
		specs = dict(nres=nres,molecule_map='one molecule many residues bigger')
		model.build_fine_to_coarse('fine','coarse_big',**specs)
	else: raise Exception('invalid molecule map %s'%molecule_map)

	# reduce each frame of theectory to coarse-grained coordinates
	nframes = len(uni.trajectory)
	frame_skip = 100
	n_cg_beads = model.reps['coarse']['length']
	coords_red = np.zeros((len(range(0,nframes,frame_skip)),n_cg_beads,3))
	for ff,fr in enumerate(range(0,nframes,frame_skip)):
		uni.trajectory[fr]
		# apply Angstroms to nm conversion here
		coords_red[ff] = model.reduce(source='fine',target='coarse',positions=sel.positions/lenscale)

	# crude method to write the coarse-grained trajectory for review
	if write_coarse_coordinates:
		"""You can review the mapping by looking at the trajectory dumped below."""
		residue_indices = model.reps['coarse']['resnums']
		atom_names = model.reps['coarse']['names']
		residue_names = np.array(['AGC' for i in residue_indices])
		for fr in range(len(coords_red)):
			polymer = GMXStructure(pts=coords_red[fr],residue_indices=residue_indices,
				residue_names=residue_names,atom_names=atom_names,box=[100.,100.,100.])
			polymer.write(state.here+'cg_model%04d.gro'%fr)
		bash('cat cg_model0* > cg_model.gro',cwd=state.here)

	#! inspecting the distributions (note that this method is used to build larger polymers)
	if inspect_distributions:
		#---! read in model here
		bead_1 = 'SB1'
		bead_2 = 'SB2'
		names = [i['name'] for i in model.reps['coarse']['beads']] 
		#---get the positions for this bond pair
		positions = coords_red[:,np.array([names.index(bead_1),names.index(bead_2)])]
		#---compute observed distances
		distances = np.linalg.norm(positions.ptp(axis=1),axis=1)
		#---! looks normal to us!
		import matplotlib as mpl;import matplotlib.pyplot as plt
		#---! check the distributions
		counts,bins = np.histogram(distances);plt.plot((bins[1:]+bins[:-1])/2.,counts);plt.show()

	# read the mapping
	with open(settings.model_spec) as fp: model_spec = yaml.load(fp.read())
	# create a model from the mapping and write dextran.itp from the MultiScaleModel class
	if molecule_map=='one molecule many residues bigger':
		model.interpret_model('coarse',coords_red,**model_spec)
		model.interpret_model('coarse_big',coords_red,**model_spec)
	elif molecule_map=='one molecule many residues injective':
		model.interpret_model('coarse',coords_red,remember_bonds=False,**model_spec)
	else: raise Exception
	#! next: run the simulation, ensure it is stable, and analyze the bond distributions

def add_antifreeze(structure,gro,ratio):
	"""Add antifreeze particles."""
	struct = GMXStructure(state.here+'%s.gro'%structure)
	inds = np.where(struct.residue_names=='W')[0]
	n_af = int(float(len(inds))*settings.antifreeze_ratio)
	subsel = np.random.permutation(range(len(inds)))[:n_af]
	struct.residue_names[inds[subsel]] = 'WF'
	struct.atom_names[inds[subsel]] = 'WF'
	component('W',count=component('W')-n_af)
	component('WF',count=n_af)
	#! requires rename_detected_composition
	struct.regroup()
	struct.write(state.here+'%s.gro'%gro)
