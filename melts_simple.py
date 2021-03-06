#!/usr/bin/env python

import os,re,subprocess
import numpy as np
import copy
import MDAnalysis
from melts_tools import vecnorm,dotplace,rotation_matrix,align,apply_rotation,review3d

_not_reported = ['vecnorm','dotplace','rotation_matrix',
	'check_cols','review3d','apply_rotation','align','vecangle','measure_torsions']

class Monomer:
	"""
	Manipulate a single monomer in a polymer. Used with lattice functions.
	"""
	termini_removals = {'start':[0],'stop':[1],'mid':[0,1]}
	def __init__(self,gro,**kwargs):
		"""
		Holds a single monomer of a growing polymer with useful metadata and placement routines.
		"""
		#---using read_gro from common
		self.raw = read_gro(gro)
		self.xyz = np.array(self.raw['points'])
		#---get the starting and ending atoms for the repeat unit vector
		self.atom_start,self.atom_end = [state.place_specs['repeat_unit'][i] for i in ['start','end']]
		#---copy the xyz points and move them to the beginning of the link
		self.monomer_origin = self.xyz[self.raw['atom_names'].index(self.atom_start)]
		self.ru_vec = self.xyz[self.raw['atom_names'].index(self.atom_end)]-self.monomer_origin
		#---rescale to the length of the repeat unit vector
		self.a0 = np.linalg.norm(self.ru_vec)
		self.remove_rules = kwargs.get('remove_rules',[])
		self.terminus = None

	def place_on_linker(self,link_vec):
		"""
		Rotate and move a monomer to a vector.
		"""
		xyzp = np.array(self.xyz) - self.monomer_origin
		#---rotate the new points so the repeat unit vector aligns with the link vector
		rotation_angle = np.arccos(np.dot(vecnorm(self.ru_vec),vecnorm(link_vec)))
		rotation_axis = np.cross(vecnorm(self.ru_vec),vecnorm(link_vec))
		rotation = rotation_matrix(rotation_axis,rotation_angle)
		xyzp_rotated = np.dot(rotation,xyzp.T).T
		return xyzp_rotated
		xyzp_placed = xyzp_rotated + polyr[mono_num] + offset
		#---if there are linkage removal rules, remove those atoms
		if removes and False:
			#---select rules for before/after atom removal
			if mono_num == 0: remove_which = [1]
			elif mono_num>0 and mono_num<len(poly)-1: remove_which = [0,1]
			else: remove_which = [0]
			for rel_i in remove_which:
				remove_inds = [monomer['atom_names'].index(r) for r in removes[rel_i]]			
				xyzp_placed = [i for ii,i in enumerate(xyzp_placed) if ii not in remove_inds]

	def lines(self,atom_num=1,resnum=1):
		"""
		Return GRO lines for this monomer. Also handles deletion rules.
		"""
		#---loop over atoms in a monomer
		lines,atom_num_abs = [],int(atom_num)
		#---figure out which deletions to apply depending on whether this monomer is a terminus
		if self.remove_rules: remove_which = self.termini_removals[self.terminus]
		else: remove_which = []
		#---select atoms that are not on the deletion list for either terminus
		atom_indices = [ii for ii,i in enumerate(self.raw['atom_names']) 
			if not any([i in self.remove_rules[j] for j in remove_which])]
		for anum in atom_indices:
			#---apply GRO formatting here
			line = (
				('%d'%resnum).rjust(5)+
				self.raw['residue_names'][anum].ljust(5)+
				self.raw['atom_names'][anum].rjust(5)+
				('%d'%atom_num_abs).rjust(5)+
				''.join([dotplace(x) for x in self.xyz[anum]])+'\n')
			lines.append(line)
			atom_num_abs += 1
		return lines

def check_cols(new_entry,column_order):
	"""Temporary function redundant with GMXTopology, used below."""
	#---! the following writer might be redundant with GMXTopology
	#---check that we have only the first contiguous sets of columns
	if len(new_entry)<len(column_order):
		if not set(new_entry.keys())==set(column_order[:len(new_entry)]):
			raise Exception('key error in columns')

def make_gel(*args,**kwargs):
	"""Switch for different kinds of gels."""
	if state.on_lattice: make_gel_on_lattice(*args,**kwargs)
	else: make_gel_off_lattice(*args,**kwargs)

def make_gel_on_lattice(name='melt',a0=1.0,sizer=10,n_p=36,volume_limit=0.2,
	uniform=True,diagonals=False,review=False):
	"""
	Make a melt by placing monomers on a (square) lattice in periodic 3-space according to a random walk.
	"""
	cwd = state.here
	#---instantiate a monomer to retrieve the correct spacing of the repeat unit (a0)
	mono = Monomer(gro='aglc.gro',remove_rules=state.place_specs.get('linkage_delete_atoms'))

	#---make a universe
	if uniform: xd,yd,zd = mono.a0,mono.a0,mono.a0
	else: raise Exception
	xld,yld,zld = sizer,sizer,sizer
	data_discrete = np.array([[i,j,k] 
		for k in np.arange(0,zld,1) 
		for j in np.arange(0,yld,1) 
		for i in np.arange(0,xld,1)]).astype(int)

	#---generate neighbor rules (with periodic boundary conditions)
	if not diagonals:
		neighbor_list = [(0+i,0+j,0+k) 
			for i in range(-1,2) for j in range(-1,2) 
			for k in range(-1,2) if not all([x==0 for x in [i,j,k]]) 
			and sum([x!=0 for x in [i,j,k]])==1]
	else: raise Exception
	neighbors = lambda x,y,z : [((x+i)%xld,(y+j)%yld,(z+k)%zld) for i,j,k in neighbor_list]

	#---CONSTRUCT A non-intersecting walk in 3-space
	#---track traffic (no intersections) and random walks through the space
	walkers,walks = [],np.zeros((xld,yld,zld)).astype(int)
	walk_id = 1
	while float(np.sum(walks==0))/np.product(walks.shape)>volume_limit:
		filled_up = np.product(walks.shape)-float(np.sum(walks==0))
		status('[CONSTRUCT] making polymers',i=filled_up,looplen=np.product(walks.shape)*(1.0-volume_limit))
		#---random starting point
		pos = tuple(np.transpose(np.where(walks==0))[np.random.randint(np.sum(walks==0))])
		#---initiate
		size = 0
		walks[pos] = walk_id
		#---loop
		this_walk = [pos]
		failure = False
		while size<n_p:
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
	status('\n[CONSTRUCT] done')

	#---use a half-lattice offset to place the dots in the center of the box if we are doing a 3D review
	if uniform: offset = a0/2.0
	else: raise Exception

	#---review in 3D and then exit (note that mayavi is asynchronous)
	if review:
		#---requires working mayavi installation (typically via qt4 backend)
		os.environ['ETS_TOOLKIT']='qt4'
		from mayavi import mlab
		import matplotlib as mpl
		import matplotlib.pylab as plt
		#---import locals
		from plotter3d import meshplot,meshpoints,pbcwire
		meshpoints(data_discrete+offset,color=(1,1,1),scale_factor=0.2,opacity=0.5)
		pbcwire([xld,yld,zld])
		for ww,walker in enumerate(walkers):
			#---identify breaks over PBCs and exclude the tube there
			breaks, = np.where(np.abs(np.sum([i[1:]-i[:-1] for i in [np.array(walker)]][0],axis=1))>1.0)
			for seg in np.split(walker,breaks+1):
				mlab.plot3d(*np.array(seg).T+offset,color=mpl.cm.jet(float(ww)/(len(walkers)-1))[:3],
					tube_radius=0.1)

	#---MOVE A MONOMER TO LINKS USING REPEAT UNIT
	combined_xyz = [[] for i in walkers]
	#---loop over all polymers
	for poly_num,poly in enumerate(walkers):
		#---loop over monomers with a polymer
		for mono_num in range(len(poly)-1):
			#---rescale the polymer (from integer units) to the correct spacing
			polyr = np.array(poly)*np.array([xd,yd,zd])
			link_vec = np.array(polyr[mono_num+1])-np.array(polyr[mono_num])
			#---the monomer can move itself into position, given a rescaled link vector
			xyz_placed = mono.place_on_linker(link_vec=link_vec)
			#---move the monomer to the start position with pbc offset
			combined_xyz[poly_num].append(xyz_placed + polyr[mono_num] + offset)

	#---view the results
	if review:
		meshpoints(xyzp_placed,scale_factor=0.1)
		meshpoints(np.concatenate((xyzp_placed[9:10],xyzp_placed[9:10])),color=(1,0,1),scale_factor=0.2)
		meshpoints(np.concatenate((xyzp_placed[22:23],xyzp_placed[22:23])),color=(1,0,1),scale_factor=0.2)
		raw_input('[QUESTION] enter anything to continue')

	#---write the GRO without atom deletions for testing purposes
	save_rules = list(mono.remove_rules)
	mono.remove_rules = []
	#---custom gro writer
	lines,resnum_abs,atom_num_abs = [],1,1
	#---loop over polymers
	for poly_num,poly in enumerate(combined_xyz):
		#---loop over monomers in each polymer
		for mono_num,mono_this in enumerate(poly):
			#---send coordinates for this monomer to the mono instance
			mono.xyz = mono_this
			#---specify the terminus and the monomer class will delete the correct molecules
			mono.terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
			lines_more = mono.lines(atom_num=atom_num_abs,resnum=resnum_abs)
			lines.extend(lines_more)
			atom_num_abs += len(lines_more)
		resnum_abs += 1
	with open(cwd+'%s_raw.gro'%name,'w') as fp: 
		fp.write('%s\n%d\n'%(name,len(lines))+''.join(lines)+' '.join(
			[dotplace(a0*i) for i in [xld,yld,zld]])+'\n')
	#---resume standard execution and repeat
	mono.remove_rules = save_rules
	#---custom gro writer
	lines,resnum_abs,atom_num_abs = [],1,1
	#---loop over polymers
	for poly_num,poly in enumerate(combined_xyz):
		#---loop over monomers in each polymer
		for mono_num,mono_this in enumerate(poly):
			#---send coordinates for this monomer to the mono instance
			mono.xyz = mono_this
			#---specify the terminus and the monomer class will delete the correct molecules
			mono.terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
			lines_more = mono.lines(atom_num=atom_num_abs,resnum=resnum_abs)
			lines.extend(lines_more)
			atom_num_abs += len(lines_more)
		resnum_abs += 1

	mol = GMXTopology(state.aglc_source)
	#---select the correct molecule from the ITP
	molspec = mol.molecules['Other']

	#---write the top file
	remove_rules = state.place_specs.get('linkage_delete_atoms')
	with open(state.here+'vacuum.top','w') as fp:
		#---refer to the CHARMM forcefield
		fp.write('#include "./charmm36.ff/forcefield.itp"\n')
		fp.write('#include "./charmm36.ff/tip3p.itp"\n')
		#---define the molecule type
		fp.write('[ moleculetype ]\n')
		fp.write('; Name nrexcl\n')
		fp.write('AGLC 3\n')
		#---write atoms
		fp.write('[ atoms ]\n')
		fp.write('; nr type resnr residue atom cgnr charge mass typeB chargeB massB\n')
		fp.write('; residue 1 AGLC rtp AGLC q  0.0\n')
		column_order = GMXTopology._entry_abstracted['atoms']['records'].split()
		#---loop over monomers in the polymer
		resnum_abs,atom_num_abs = 1,1
		for mono_num in range(n_p):
			#---loop over atoms in a monomer
			#---figure out which deletions to apply depending on whether this monomer is a terminus
			terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
			remove_which = Monomer.termini_removals[terminus]
			#---select atoms that are not on the deletion list for either terminus
			atom_indices = [ii for ii,i in enumerate(molspec['atoms']) 
				if not any([i['atom'] in remove_rules[j] for j in remove_which])]
			#---write atom entries for the valid atom_indices
			for atom_num in atom_indices:
				new_entry = dict(molspec['atoms'][atom_num])
				#---update atom number and residue index
				new_entry['id'] = atom_num_abs
				new_entry['resnr'] = mono_num+1
				new_entry['cgnr'] = atom_num_abs
				#---loop over possible atom name/partial charge changes
				if state.place_specs['atom_name_changes']:
					assert type(state.place_specs['atom_name_changes'])==list,'changes must be a list'
					for change_entry in state.place_specs['atom_name_changes']:
						#---rule is a "next" and the monomer is not the first one
						#---...or rule is a "previous" and the monomer is not the last one
						if (((mono_num>0 and change_entry['which'] == 'next') or 
							(mono_num<n_p-1 and change_entry['which'] == 'previous'))
							and change_entry['atom']==new_entry['atom']):
							#---rule is a "next" and the monomer is not the first one
							#---...or rule is a "previous" and the monomer is not the last one
							assert new_entry['type']==change_entry['from']
							new_entry['type'] = change_entry['to']
							#---check the current partial charges to prevent user error
							if not float(new_entry['charge'])==float(change_entry['from_charge']):
								raise Exception('previous (partial) charge does not equal the new charge '+
									'please check the values and/or modify the precision to pass this test')
							#---fix the partial charges
							new_entry['charge'] = '%.3f'%float(change_entry['to_charge'])
				#---use the column order to print the entry correctly
				check_cols(new_entry,column_order)
				line = ' '.join([str(new_entry[i]) for i in column_order[:len(new_entry)]])
				fp.write(line+'\n')
				atom_num_abs += 1

		#---write bond entries
		fp.write('[ bonds ]\n')
		fp.write('; ai aj funct c0 c1 c2 c3\n')
		column_order = GMXTopology._entry_abstracted['bonds']['records'].split()
		#---loop over monomers in the polymer
		previous_atoms,linker_indices,all_bonds = 0,[],[]
		for mono_num in range(n_p):
			#---! the following terminus code is redundant -- make it a function
			#---figure out which deletions to apply depending on whether this monomer is a terminus
			terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
			remove_which = Monomer.termini_removals[terminus]
			#---select atoms that are not on the deletion list for either terminus
			atom_indices = [ii for ii,i in enumerate(molspec['atoms']) 
				if not any([i['atom'] in remove_rules[j] for j in remove_which])]
			#---add the linker if necessary
			if mono_num > 0:
				these_atoms = [i['atom'] for i in molspec['atoms']]
				#---get the index of the "end" of the repeat unit
				prev_index_abs = these_atoms.index(state.place_specs['repeat_unit']['end'])
				link_start = atom_indices_prev.index(prev_index_abs) + previous_atoms_prev + 1
				link_end = [these_atoms[j] for j in atom_indices].index(
					state.place_specs['repeat_unit']['next']) + 1 + previous_atoms
				#---! manually specifying function here
				linker_bond = {'i':'%s'%link_start,'length':'','j':'%s'%link_end,'force': '1','funct':''}
				linker_indices.append((link_start,link_end))
				line = ' '.join([str(linker_bond[i]) for i in column_order])
				fp.write(line+'\n')
			#---only include bonds for atoms that are not removed
			valid_bonds = [i for i in molspec['bonds'] if not any([int(i[j])-1 
				not in atom_indices for j in 'ij'])]
			#---loop over the bonds and increment residue indices
			for bond_num,entry in enumerate(valid_bonds):
				new_entry = dict(entry)
				#---update indices
				new_entry['i'] = atom_indices.index(int(new_entry['i'])-1)+1 + previous_atoms
				new_entry['j'] = atom_indices.index(int(new_entry['j'])-1)+1 + previous_atoms
				all_bonds.append((new_entry['i'],new_entry['j']))
				#---use the column order to print the entry correctly
				check_cols(new_entry,column_order)
				line = ' '.join([str(new_entry[i]) for i in column_order[:len(new_entry)]])
				fp.write(line+'\n')
			previous_atoms_prev = int(previous_atoms)
			previous_atoms += len(atom_indices)
			#---retain these atom indices for the linker
			atom_indices_prev = list(atom_indices)

		#---write ANGLE entries
		fp.write('[ angles ]\n')
		fp.write('; i j k funct angle force\n')
		column_order = GMXTopology._entry_abstracted['angles']['records'].split()
		#---loop over monomers in the polymer
		previous_atoms = 0
		for mono_num in range(n_p):
			#---! the following terminus code is redundant -- make it a function
			#---figure out which deletions to apply depending on whether this monomer is a terminus
			terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
			#terminus = {0:'stop',n_p-1:'start'}.get(mono_num,'mid')
			remove_which = Monomer.termini_removals[terminus]
			#---select atoms that are not on the deletion list for either terminus
			atom_indices = [ii for ii,i in enumerate(molspec['atoms']) 
				if not any([i['atom'] in remove_rules[j] for j in remove_which])]
			#---only include bonds for atoms that are not removed
			valid_angles = [i for i in molspec['angles'] if not any([int(i[j])-1 
				not in atom_indices for j in 'ijk'])]
			#---loop over the bonds and increment residue indices
			for angle_num,entry in enumerate(valid_angles):
				new_entry = dict(entry)
				#---update indices
				new_entry['i'] = atom_indices.index(int(new_entry['i'])-1)+1 + previous_atoms
				new_entry['j'] = atom_indices.index(int(new_entry['j'])-1)+1 + previous_atoms
				new_entry['k'] = atom_indices.index(int(new_entry['k'])-1)+1 + previous_atoms
				#---use the column order to print the entry correctly
				check_cols(new_entry,column_order)
				line = ' '.join([str(new_entry[i]) for i in column_order[:len(new_entry)]])
				fp.write(line+'\n')
			previous_atoms_prev = int(previous_atoms)
			previous_atoms += len(atom_indices)
			#---retain these atom indices for the linker
			atom_indices_prev = list(atom_indices)
			#---add angles for the linker
			for lnum,(linkl,linkr) in enumerate(linker_indices):
				bound_not = lambda p : [j for k in [i for i in all_bonds if p in i] for j in k if j!=p]
				l1,r1 = [list(set(bound_not(p))) for p in [linkl,linkr]]
				new_angles = [(linkl,linkr,i) for i in r1]+[(i,linkl,linkr) for i in l1]
				for na in new_angles:
					new_entry = dict(entry)
					for ll,l in enumerate('ijk'): new_entry[l] = na[ll]
					#---use the column order to print the entry correctly
					check_cols(new_entry,column_order)
					line = ' '.join([str(new_entry[i]) for i in column_order[:len(new_entry)]])
					fp.write(line+'\n')

		#---write DIHEDRAL entries
		fp.write('[ dihedrals ]\n')
		fp.write('; i j k l funct angle force\n')
		column_order = GMXTopology._entry_abstracted['dihedrals']['records'].split()
		#---loop over monomers in the polymer
		previous_atoms = 0
		for mono_num in range(n_p):
			#---! the following terminus code is redundant -- make it a function
			#---figure out which deletions to apply depending on whether this monomer is a terminus
			terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
			remove_which = Monomer.termini_removals[terminus]
			#---select atoms that are not on the deletion list for either terminus
			atom_indices = [ii for ii,i in enumerate(molspec['atoms']) 
				if not any([i['atom'] in remove_rules[j] for j in remove_which])]
			#---only include bonds for atoms that are not removed
			valid_dihedrals = [i for i in molspec['dihedrals'] if not any([int(i[j])-1 
				not in atom_indices for j in 'ijkl'])]
			#---loop over the bonds and increment residue indices
			for angle_num,entry in enumerate(valid_dihedrals):
				new_entry = dict(entry)
				#---update indices
				new_entry['i'] = atom_indices.index(int(new_entry['i'])-1)+1 + previous_atoms
				new_entry['j'] = atom_indices.index(int(new_entry['j'])-1)+1 + previous_atoms
				new_entry['k'] = atom_indices.index(int(new_entry['k'])-1)+1 + previous_atoms
				new_entry['l'] = atom_indices.index(int(new_entry['l'])-1)+1 + previous_atoms
				#---use the column order to print the entry correctly
				check_cols(new_entry,column_order)
				line = ' '.join([str(new_entry[i]) for i in column_order[:len(new_entry)]])
				fp.write(line+'\n')
			previous_atoms_prev = int(previous_atoms)
			previous_atoms += len(atom_indices)
			#---retain these atom indices for the linker
			atom_indices_prev = list(atom_indices)
		#---add dihedrals for the linker
		for lnum,(linkl,linkr) in enumerate(linker_indices):
			bound_not = lambda p : [j for k in [i for i in all_bonds if p in i] for j in k if j!=p]
			l1,r1 = [list(set(bound_not(p))) for p in [linkl,linkr]]
			l2 = [(i,j,linkl,linkr) for j in l1 for i in bound_not(j) if i!=linkl]
			r2 = [(linkl,linkr,j,i) for j in r1 for i in bound_not(j) if i!=linkr]
			#---still need the ones that run all the way through
			thru = [(l,linkl,linkr,r) for l in l1 for r in r1]
			combos = r2+l2+thru
			#---no reversed dihedrals
			extra_dihedrals = [combos[k] for k in [ii for ii,i in enumerate(combos) 
				if i not in [list(j)[::-1] for j in combos]]]
			for ed in extra_dihedrals:
				new_entry = dict(entry)
				for ll,l in enumerate('ijkl'): new_entry[l] = ed[ll]
				#---use the column order to print the entry correctly
				check_cols(new_entry,column_order)
				line = ' '.join([str(new_entry[i]) for i in column_order[:len(new_entry)]])
				fp.write(line+'\n')
		#---write the total number of molecules
		fp.write('[ system ]\nMY SYSTEM\n[ molecules ]\nAGLC %d\n'%len(combined_xyz))

	#---write the GRO
	with open(cwd+'%s.gro'%name,'w') as fp: 
		fp.write('%s\n%d\n'%(name,len(lines))+''.join(lines)+' '.join(
			[dotplace(a0*i) for i in [xld,yld,zld]])+'\n')

	state.n_polymers = len(combined_xyz)
	#---! commands that trail this function in the dextran development code
	for i in state.include_adds: include(i)
	#---! note that the QTOT may be incorrect here. check it!
	component(state.polymer_name,count=state.n_polymers)

def make_off_lattice_walk(n_segments,length,angle,torsion):
	"""
	Make an unconstrained 3D random walk given a link length a0 and a mean angle and dihedral.
	Angle in radians. Note that "dihedral" really means torsion, since this might be a coarse-model.
	"""
	angle_rad = angle/180*np.pi
	torsion_rad = torsion/180*np.pi
	#---! enforce ranges on the angles and dihedral
	if n_segments<3: raise Exception('this walk requires an angle and a dihedral (min. 3 links)')
	pts = np.zeros((n_segments+1,3))
	for mnum in range(n_segments+1):
		if mnum==0: pts[mnum] = [0,0,0]
		elif mnum==1: pts[mnum] = [length,0,0]
		elif mnum==2: 
			rotation = rotation_matrix([0.,0.,-1.],angle_rad)
			rotated_pt = np.dot(rotation,np.array([pts[mnum-2]-pts[mnum-1]]).T).T
			pts[mnum] = vecnorm(rotated_pt+pts[mnum-1])*length
		elif mnum>=3: 
			#---given pts 1,2,3 and angle,dihedral then find point 4
			#---get the vector normal to vectors 1-2,2-3
			vec_12x23 = vecnorm(np.cross(vecnorm(pts[mnum-3]-pts[mnum-2]),vecnorm(pts[mnum-1]-pts[mnum-2])))
			vec_23 = vecnorm(pts[mnum-2]-pts[mnum-1])
			#---pts4a is normal to 23 and to the vector normal to 12 and 23
			pts4a = np.cross(vec_12x23,vec_23) + pts[mnum-1]
			#---now we wish to find pts4b, which is 4a rotated in the plane normal to 23 so that it has the
			#---...correct dihedral
			#---we reverse the sign of the dihedral to obey a right-hand rule
			#---...this means that if you are looking down the 23 vector of the dihedral, a positive 
			#---...dihedral angle causes the 4th point to rotate clockwise
			rotation = rotation_matrix(vec_23,-1*torsion_rad)
			pts4b = np.dot(rotation,np.array([pts4a-pts[mnum-1]]).T).T[0] + pts[mnum-1]
			#---rotate pts4b in the plane normal to the vector that is normal to the vector from the third
			#---...point in the dihedral (pts[mnum]-1) to pts4b (which has the right dihedral angle) and also
			#---...normal to the 23 vector. this ensures that the current rotation does not change the 
			#---...torsion. we first subtract 90 so that the angle is applied relative to the 23 vector
			vec_34b = np.cross(vec_23,vecnorm(pts4b-pts[mnum-1]))
			rotation = rotation_matrix(vec_34b,angle_rad-np.pi/2.)
			pts4c = vecnorm(np.dot(rotation,np.array([pts4b-pts[mnum-1]]).T).T[0])*length + pts[mnum-1]
			pts[mnum] = pts4c
	if state.review3d: review3d(lines=[pts],points=[pts],tube=0.02,radius=0.2)
	return pts

def make_gel_off_lattice(name='melt',a0=1.0,sizer=10,n_p=36,volume_limit=0.2,
	uniform=True,diagonals=False,review=False,cwd=None,angle=90.,torsion=90.):
	"""
	Make a melt by placing monomers on a (square) lattice in periodic 3-space according to a random walk.
	"""
	cwd = state.here
	n_p = state.melt_settings['n_p']

	#---instantiate a monomer to retrieve the correct spacing of the repeat unit (a0)
	mono = Monomer(gro='aglc.gro',remove_rules=state.place_specs.get('linkage_delete_atoms'))

	#---prepare an abstract 3D walk
	walk_abstract_pts = make_off_lattice_walk(
		n_p,state.melt_settings['a0'],
		angle,torsion)

	#---loop over "links" in our 3D walk and superimpose a monomer on each link
	points_raw = []
	origin = np.array([0,0,0])
	for mnum in range(1,n_p+1):
		link_vec = walk_abstract_pts[mnum]-walk_abstract_pts[mnum-1]
		xyz_placed = mono.place_on_linker(link_vec=link_vec)
		points_raw.append(xyz_placed + walk_abstract_pts[mnum-1])
	points_raw = np.array(points_raw)
	if state.review3d: review3d(points=[points_raw])

	#---! temporarily back to "combined_xyz" nomenclature
	combined_xyz = [points_raw]

	#---! HACKED
	box_vectors = [100.,100.,100.]

	#---write the GRO without atom deletions for testing purposes
	save_rules = list(mono.remove_rules)
	mono.remove_rules = []
	#---custom gro writer
	lines,resnum_abs,atom_num_abs = [],1,1
	#---loop over polymers
	for poly_num,poly in enumerate(combined_xyz):
		#---loop over monomers in each polymer
		for mono_num,mono_this in enumerate(poly):
			#---send coordinates for this monomer to the mono instance
			mono.xyz = mono_this
			#---specify the terminus and the monomer class will delete the correct molecules
			mono.terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
			lines_more = mono.lines(atom_num=atom_num_abs,resnum=resnum_abs)
			lines.extend(lines_more)
			atom_num_abs += len(lines_more)
		resnum_abs += 1
	#---! eventually remove the raw form because we don't need to see that version
	with open(cwd+'%s_raw.gro'%name,'w') as fp: 
		fp.write('%s\n%d\n'%(name,len(lines))+''.join(lines)+' '.join(
			[dotplace(a0*i) for i in box_vectors])+'\n')
	#---resume standard execution and repeat
	mono.remove_rules = save_rules
	#---custom gro writer
	lines,resnum_abs,atom_num_abs = [],1,1
	#---loop over polymers
	for poly_num,poly in enumerate(combined_xyz):
		#---loop over monomers in each polymer
		for mono_num,mono_this in enumerate(poly):
			#---send coordinates for this monomer to the mono instance
			mono.xyz = mono_this
			#---specify the terminus and the monomer class will delete the correct molecules
			mono.terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
			lines_more = mono.lines(atom_num=atom_num_abs,resnum=resnum_abs)
			lines.extend(lines_more)
			atom_num_abs += len(lines_more)
		resnum_abs += 1
	#---write the GRO
	with open(cwd+'%s.gro'%name,'w') as fp: 
		fp.write('%s\n%d\n'%(name,len(lines))+''.join(lines)+' '.join(
			[dotplace(a0*i) for i in box_vectors])+'\n')

	#---! example topology via: topex = GMXTopology('inputs/test.top')
	#---note that the on-lattice version had a custom topology writer. here we try the automatic way.
	top = GMXTopology()

	#---construct the new molecule
	new_mol = {'moleculetype':{'molname':'AGLC','nrexcl':3},
		'atoms':[],'bonds':[],'angles':[],'dihedrals':[]}

	mol = GMXTopology(state.aglc_source)
	#---select the correct molecule from the ITP
	#---! again note that the QTOT may be incorrect and needs checked
	molspec = mol.molecules['Other']

	#---build atoms
	resnum_abs,atom_num_abs = 1,1
	for mono_num in range(n_p):
		#---loop over atoms in a monomer
		#---figure out which deletions to apply depending on whether this monomer is a terminus
		terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
		remove_which = Monomer.termini_removals[terminus]
		#---select atoms that are not on the deletion list for either terminus
		atom_indices = [ii for ii,i in enumerate(molspec['atoms']) 
			if not any([i['atom'] in mono.remove_rules[j] for j in remove_which])]
		#---write atom entries for the valid atom_indices
		for atom_num in atom_indices:
			new_entry = dict(molspec['atoms'][atom_num])
			#---update atom number and residue index
			new_entry['id'] = atom_num_abs
			new_entry['resnr'] = mono_num+1
			new_entry['cgnr'] = atom_num_abs
			#---loop over possible atom name/partial charge changes
			#---! isn't the following always true? redundant?
			if state.place_specs['atom_name_changes']:
				assert type(state.place_specs['atom_name_changes'])==list,'changes must be a list'
				for change_entry in state.place_specs['atom_name_changes']:
					#---rule is a "next" and the monomer is not the first one
					#---...or rule is a "previous" and the monomer is not the last one
					if (((mono_num>0 and change_entry['which'] == 'next') or 
						(mono_num<n_p-1 and change_entry['which'] == 'previous'))
						and change_entry['atom']==new_entry['atom']):
						#---rule is a "next" and the monomer is not the first one
						#---...or rule is a "previous" and the monomer is not the last one
						assert new_entry['type']==change_entry['from']
						new_entry['type'] = change_entry['to']
						#---check the current partial charges to prevent user error
						if not float(new_entry['charge'])==float(change_entry['from_charge']):
							raise Exception('previous (partial) charge does not equal the new charge '+
								'please check the values and/or modify the precision to pass this test')
						#---fix the partial charges
						new_entry['charge'] = '%.3f'%float(change_entry['to_charge'])
			#---add the new entry to our new molecule
			new_mol['atoms'].append(new_entry)
			atom_num_abs += 1

	#---build bonds
	previous_atoms,linker_indices,all_bonds = 0,[],[]
	for mono_num in range(n_p):
		#---! the following terminus code is redundant -- make it a function
		#---figure out which deletions to apply depending on whether this monomer is a terminus
		terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
		remove_which = Monomer.termini_removals[terminus]
		#---select atoms that are not on the deletion list for either terminus
		atom_indices = [ii for ii,i in enumerate(molspec['atoms']) 
			if not any([i['atom'] in mono.remove_rules[j] for j in remove_which])]
		#---add the linker if necessary
		if mono_num > 0:
			these_atoms = [i['atom'] for i in molspec['atoms']]
			#---get the index of the "end" of the repeat unit
			prev_index_abs = these_atoms.index(state.place_specs['repeat_unit']['end'])
			link_start = atom_indices_prev.index(prev_index_abs) + previous_atoms_prev + 1
			link_end = [these_atoms[j] for j in atom_indices].index(
				state.place_specs['repeat_unit']['next']) + 1 + previous_atoms
			#---! manually specifying function here
			linker_bond = {'i':'%s'%link_start,'length':'','j':'%s'%link_end,'force': '1','funct':''}
			linker_indices.append((link_start,link_end))
			new_mol['bonds'].append(linker_bond)
		#---only include bonds for atoms that are not removed
		valid_bonds = [i for i in molspec['bonds'] if not any([int(i[j])-1 
			not in atom_indices for j in 'ij'])]
		#---loop over the bonds and increment residue indices
		for bond_num,entry in enumerate(valid_bonds):
			new_entry = dict(entry)
			#---update indices
			new_entry['i'] = atom_indices.index(int(new_entry['i'])-1)+1 + previous_atoms
			new_entry['j'] = atom_indices.index(int(new_entry['j'])-1)+1 + previous_atoms
			all_bonds.append((new_entry['i'],new_entry['j']))
			new_mol['bonds'].append(new_entry)
		previous_atoms_prev = int(previous_atoms)
		previous_atoms += len(atom_indices)
		#---retain these atom indices for the linker
		atom_indices_prev = list(atom_indices)

	#---build angles
	previous_atoms = 0
	for mono_num in range(n_p):
		#---! the following terminus code is redundant -- make it a function
		#---figure out which deletions to apply depending on whether this monomer is a terminus
		terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
		#terminus = {0:'stop',n_p-1:'start'}.get(mono_num,'mid')
		remove_which = Monomer.termini_removals[terminus]
		#---select atoms that are not on the deletion list for either terminus
		atom_indices = [ii for ii,i in enumerate(molspec['atoms']) 
			if not any([i['atom'] in mono.remove_rules[j] for j in remove_which])]
		#---only include bonds for atoms that are not removed
		valid_angles = [i for i in molspec['angles'] if not any([int(i[j])-1 
			not in atom_indices for j in 'ijk'])]
		#---loop over the bonds and increment residue indices
		for angle_num,entry in enumerate(valid_angles):
			new_entry = dict(entry)
			#---update indices
			new_entry['i'] = atom_indices.index(int(new_entry['i'])-1)+1 + previous_atoms
			new_entry['j'] = atom_indices.index(int(new_entry['j'])-1)+1 + previous_atoms
			new_entry['k'] = atom_indices.index(int(new_entry['k'])-1)+1 + previous_atoms
			new_mol['angles'].append(new_entry)
		previous_atoms_prev = int(previous_atoms)
		previous_atoms += len(atom_indices)
		#---retain these atom indices for the linker
		atom_indices_prev = list(atom_indices)
		#---add angles for the linker
		for lnum,(linkl,linkr) in enumerate(linker_indices):
			bound_not = lambda p : [j for k in [i for i in all_bonds if p in i] for j in k if j!=p]
			l1,r1 = [list(set(bound_not(p))) for p in [linkl,linkr]]
			new_angles = [(linkl,linkr,i) for i in r1]+[(i,linkl,linkr) for i in l1]
			for na in new_angles:
				new_entry = dict(entry)
				for ll,l in enumerate('ijk'): 
					new_entry[l] = na[ll]
					new_mol['angles'].append(new_entry)

	#---build dihedrals
	previous_atoms = 0
	for mono_num in range(n_p):
		#---! the following terminus code is redundant -- make it a function
		#---figure out which deletions to apply depending on whether this monomer is a terminus
		terminus = {0:'start',n_p-1:'stop'}.get(mono_num,'mid')
		remove_which = Monomer.termini_removals[terminus]
		#---select atoms that are not on the deletion list for either terminus
		atom_indices = [ii for ii,i in enumerate(molspec['atoms']) 
			if not any([i['atom'] in mono.remove_rules[j] for j in remove_which])]
		#---only include bonds for atoms that are not removed
		valid_dihedrals = [i for i in molspec['dihedrals'] if not any([int(i[j])-1 
			not in atom_indices for j in 'ijkl'])]
		#---loop over the bonds and increment residue indices
		for angle_num,entry in enumerate(valid_dihedrals):
			new_entry = dict(entry)
			#---update indices
			new_entry['i'] = atom_indices.index(int(new_entry['i'])-1)+1 + previous_atoms
			new_entry['j'] = atom_indices.index(int(new_entry['j'])-1)+1 + previous_atoms
			new_entry['k'] = atom_indices.index(int(new_entry['k'])-1)+1 + previous_atoms
			new_entry['l'] = atom_indices.index(int(new_entry['l'])-1)+1 + previous_atoms
			new_mol['dihedrals'].append(new_entry)
		previous_atoms_prev = int(previous_atoms)
		previous_atoms += len(atom_indices)
		#---retain these atom indices for the linker
		atom_indices_prev = list(atom_indices)
	#---add dihedrals for the linker
	for lnum,(linkl,linkr) in enumerate(linker_indices):
		bound_not = lambda p : [j for k in [i for i in all_bonds if p in i] for j in k if j!=p]
		l1,r1 = [list(set(bound_not(p))) for p in [linkl,linkr]]
		l2 = [(i,j,linkl,linkr) for j in l1 for i in bound_not(j) if i!=linkl]
		r2 = [(linkl,linkr,j,i) for j in r1 for i in bound_not(j) if i!=linkr]
		#---still need the ones that run all the way through
		thru = [(l,linkl,linkr,r) for l in l1 for r in r1]
		combos = r2+l2+thru
		#---no reversed dihedrals
		extra_dihedrals = [combos[k] for k in [ii for ii,i in enumerate(combos) 
			if i not in [list(j)[::-1] for j in combos]]]
		for ed in extra_dihedrals:
			new_entry = dict(entry)
			for ll,l in enumerate('ijkl'): 
				new_entry[l] = ed[ll]
				new_mol['dihedrals'].append(new_entry)

	#---add the molecule and write the itp
	top.add_molecule(**{'AGLC':new_mol})
	top.write(state.here+'dextran.itp')
	if not state.itp: state.itp = []
	state.itp.append('dextran.itp')

	#---collect and save the metadata
	####state.include_adds = ['./charmm36.ff/forcefield.itp','./charmm36.ff/tip3p.itp','dextran.itp']
	state.n_polymers = len(combined_xyz)
	###for i in state.include_adds: include(i)
	component(state.polymer_name,count=state.n_polymers)
	
	write_top('vacuum.top')

def make_cg_gel(name='melt',**kwargs):
	"""
	Make a melt by placing monomers on a (square) lattice in periodic 3-space according to a random walk.
	"""
	polymer_molame = 'DEX'
	monomer_resname = 'AGLC'
	cwd = state.here
	n_p = state.melt_settings['n_p']
	angle,torsion = state.melt_settings['angle'],state.melt_settings['torsion']
	a0 = state.melt_settings['a0']
	
	#---prepare an abstract 3D walk
	walk_abstract_pts = make_off_lattice_walk(
		n_p,a0,angle,torsion)
	points_with_sides = []
	for pt in walk_abstract_pts:
		points_with_sides.append(pt)
		for i in range(2):
			#---! randomly place sidechain points 1 a0 distance away. map 0 to 1 to -1 to 1
			points_with_sides.append(pt+(1-2*np.random.rand(3))/10.)
	points_with_sides = np.array(points_with_sides)
	residue_names = np.array([monomer_resname for p in points_with_sides])
	residue_indices = np.array([i/3+1 for i in range((n_p+1)*3)])
	atom_names = np.tile(['B1','B2','B3'],n_p+1)
	polymer = GMXStructure(pts=points_with_sides,residue_indices=residue_indices,
		residue_names=residue_names,atom_names=atom_names,box=[10.,10.,10.])
	polymer.write(state.here+'%s-built.gro'%name)
	gmx('editconf',structure='%s-built'%name,gro=name,c=True,log='editconf-center-polymer')

	top = GMXTopology()
	new_mol = {'bonds':[],'angles':[],'dihedrals':[],'moleculetype':{'molname':polymer_molame,'nrexcl':3}}
	#---monomer spec
	monomer = [
		{'type':'P1','mass':'60.0528','atom':'B1'},
		{'type':'P4','mass':'60.0528','atom':'B2'},
		{'type':'P4','mass':'60.0528','atom':'B3'}]
	#---write atoms
	new_mol['atoms'] = []
	atom_counter = 1
	for mnum in range(n_p+1):
		new_monomer = copy.deepcopy(monomer)
		for anum,atom in enumerate(new_monomer):
			new_monomer[anum]['id'] = atom_counter
			atom_counter += 1
			new_monomer[anum]['resnr'] = new_monomer[anum]['cgnr'] = mnum+1
			new_monomer[anum]['charge'] = '0'
			new_monomer[anum]['resname'] = monomer_resname
		new_mol['atoms'].extend(new_monomer)

	#---SECTION 1
	#---! extremely crude way to write bonds -- systematize this later!
	for mnum in range(n_p+1-1):
		#---1-4 bond
		new_mol['bonds'].append({'i':mnum*3+1,'length':'0.598','j':(mnum+1)*3+1,'force':'3000','funct':'1'})
	
	#---SECTION 2
	for mnum in range(n_p+1):
		#---1-2 and 1-3 bonds
		new_mol['bonds'].append({'i':mnum*3+1,'length':'0.598','j':mnum*3+1+1,'force':'3000','funct':'1'})
		new_mol['bonds'].append({'i':mnum*3+1,'length':'0.598','j':mnum*3+1+2,'force':'3000','funct':'1'})

	#---SECTION 3
	#---! extremely crude method for angles
	angle_template = {'i':0,'j':0,'k':0,'funct':2,'angle':180.0,'force':2000.0}
	for mnum in range(n_p+1-2):
		new_mol['angles'].append(dict(angle_template,i=mnum*3+1,j=(mnum+1)*3+1,k=(mnum+2)*3+1))

	#---SECTION 4
	#---! just a placeholder: 1-2-3 angle is 180 -- this keeps the sidechains in a straight line across back
	angle_template_flank = {'i':0,'j':0,'k':0,'funct':2,'angle':180.0,'force':2000.0}
	for mnum in range(n_p+1):
		new_mol['angles'].append(dict(angle_template_flank,i=mnum*3+1+1,j=mnum*3+1+0,k=mnum*3+1+2))

	#---SECTION 5
	#---! just a placeholder: 1-2-4 angle which orients the side chain normal to the backbone
	angle_template_flank_ortho = {'i':0,'j':0,'k':0,'funct':2,'angle':90.0,'force':2000.0}
	for mnum in range(n_p+1-1):
		new_mol['angles'].append(dict(angle_template_flank_ortho,i=mnum*3+1+1,j=mnum*3+1+0,k=(mnum+1)*3+1))
		new_mol['angles'].append(dict(angle_template_flank_ortho,i=mnum*3+1+2,j=mnum*3+1+0,k=(mnum+1)*3+1))

	#---SECTION 6
	if False:
		#---! extremely crude method for dihedrals
		dihedral_template = {'i':0,'j':0,'k':0,'l':0,'funct':1,'angle':142.0,'force':2000.0,'multiplicity':1}
		for mnum in range(n_p+1-3):
			new_mol['dihedrals'].append(dict(dihedral_template,
				i=mnum*3+1,j=(mnum+1)*3+1,k=(mnum+2)*3+1,l=(mnum+3)*3+1))

	"""
	history of manipulating which bonds get included (2017.06.15):
		started with a previous run that only had sections 1-4
			and the backbone was pretty straight 
			but the sidechains were wobbly (see "arch-v024-before-adding-124-angles")
		while reviewing the method, we added section 5 
			and noticed that the sidechains were less wobbly
			for some reason the dihedrals of the sidechains looked really straight
			but there were no dihedrals in the itp file. weird.
		added section 6 and made the torsion for the backbone 142.0 from the melts_tuner.py
	"""

	top.add_molecule(**{polymer_molame:new_mol})
	top.write(state.here+'dextran.itp')
	state.itp = ['dextran.itp']
	component(polymer_molame,count=1)
