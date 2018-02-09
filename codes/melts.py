#!/usr/bin/env python

import os,re,subprocess
import numpy as np
import copy
import MDAnalysis

_not_reported = ['vecnorm','dotplace','rotation_matrix','check_cols','review3d','apply_rotation','align']
#---! somehow status was getting reported when I merged in melts_tuner.py
_not_reported += ['status','vecangle','measure_torsions']

###---TOOLS

vecnorm = lambda x: x/np.linalg.norm(x)
dotplace = lambda n: re.compile(r'(\d)0+$').sub(r'\1',"%8.3f"%float(n)).ljust(8)

def rotation_matrix(axis,theta):
	"""
	Return the rotation matrix associated with counterclockwise rotation about
	the given axis by theta radians using Euler-Rodrigues formula.
	"""
	axis = np.asarray(axis)
	theta = np.asarray(theta)
	if all(axis==0): return np.identity(3) 
	axis = axis/np.sqrt(np.dot(axis,axis))
	a = np.cos(theta/2)
	b,c,d = -axis*np.sin(theta/2)
	aa,bb,cc,dd = a*a,b*b,c*c,d*d
	bc,ad,ac,ab,bd,cd = b*c,a*d,a*c,a*b,b*d,c*d
	return np.array([[aa+bb-cc-dd,2*(bc+ad),2*(bd-ac)],[2*(bc-ad),aa+cc-bb-dd,2*(cd+ab)],
		[2*(bd+ac),2*(cd-ab),aa+dd-bb-cc]])

def align(ref,target):
	"""Perform an alignment relative to the reference (first) structure."""
	r0,r1 = np.array(ref),np.array(target)
	origin = r0.mean(axis=0)
	r0 -= origin
	r1 -= r1.mean(axis=0)
	U,s,Vt = np.linalg.svd(np.dot(r0.T,r1))
	signer = np.identity(3)
	signer[2,2] = np.sign(np.linalg.det(np.dot(Vt.T,U)))
	RM = np.dot(np.dot(U,signer),Vt)
	r1p = apply_rotation(r1,RM)
	rmsd = np.sqrt(np.mean(np.sum((r0.T-np.dot(RM,r1.T))**2,axis=0)))
	return {'new':r1p,'rm':RM,'origin':origin,'rmsd':rmsd}

def apply_rotation(pts,matrix):
	"""Standard way to apply a rotation (possibly calculated from align)."""
	return np.dot(matrix,pts.T).T

###---VISUALIZATION

def review3d(**kwargs):
	"""
	Review things in 3D. Originally developed for coiled coils.
	!Consolidate from structural_biology.py to a central location with a standard interface. 
	Note that scaling factors assume nanometers.
	"""
	import mayavi
	from mayavi import mlab
	def pointplot(x,y,z,colormap='Spectral',r=1.0,color=(0,0,0)):
		mlab.points3d(x,y,z,colormap='Spectral',scale_factor=r,color=color)
	def lineplot(x,y,z,colormap='Spectral',r=0.1,color=(0,0,0)):
		mlab.plot3d(x,y,z,tube_radius=r,color=color)
	tube_thickness = kwargs.get('tube',0.1)
	sphere_radius = kwargs.get('radius',0.1)
	for lines in kwargs.get('lines',[]): lineplot(*lines.T,r=tube_thickness)
	for points in kwargs.get('points',[]): pointplot(*points.T,r=sphere_radius)
	mlab.show()

###---MONOMER CLASS

class Monomer:
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
		#---! HOW ARE PBCS IMPLEMENTED HERE?
		#if not np.abs(np.linalg.norm(link_vec)-self.a0)<=0.001:
		#	return 'aaa'
		#	raise Exception('[ERROR] the link vector must be the same length as the monomer repeat unit')
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
	"""
	Temporary function redundant with GMXTopology
	!replace this if you reformulate the topology stuff below.
	!remove from _not_reported when you remove it
	"""
	#---! the following writer might be redundant with GMXTopology
	#---check that we have only the first contiguous sets of columns
	if len(new_entry)<len(column_order):
		if not set(new_entry.keys())==set(column_order[:len(new_entry)]):
			raise Exception('key error in columns')

###---BASIC GEL CONSTRUCTION

def make_gel(*args,**kwargs):
	if state.on_lattice: make_gel_on_lattice(*args,**kwargs)
	else: make_gel_off_lattice(*args,**kwargs)

def make_gel_on_lattice(name='melt',a0=1.0,sizer=10,n_p=36,volume_limit=0.2,
	uniform=True,diagonals=False,review=False):

	"""
	Make a melt by placing monomers on a (square) lattice in periodic 3-space according to a random walk.
	"""

	#---destination
	#cwd = './' if not cwd else os.path.join(cwd,'')
	cwd = state.here

	#---instantiate a monomer to retrieve the correct spacing of the repeat unit (a0)
	mono = Monomer(gro='aglc.gro',remove_rules=state.place_specs.get('linkage_delete_atoms'))
	#mono = Monomer(gro='aglc.gro')

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

	#---! ongoing off_lattice development should start here

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

	#---PROTOTYPE POLYMER TOPOLOGY BUILDER
	#mol = MoleculeModel('inputs/dev-melt-aglc-charmm/aglc-example.top')
	mol = GMXTopology(state.aglc_source)
	#---select the correct molecule from the ITP
	molspec = mol.molecules['Other']
	#---! the above needs to be automated, change "Other" to AGLC 
	#---! ...and use pdb2gmx to auto-generate the monomer topology

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
	component(state.polymer_name,count=state.n_polymers)
	#---!!! THE QTOT IS WRONG AND WHEN YOU FIX THE TOP FILE YOU SHOULD FIX IT!!!!!

def make_off_lattice_walk(n_segments,length,angle,torsion):
	"""
	Make an unconstrained 3D random walk given a link length a0 and a mean angle and dihedral.
	Angle in radians.

	psuedocode:
		start from the origin
		first linker is along x-axis, 1 unit lenth
		second linker is the incoming angle, rotated in e.g. xy-plane
		third linker must have a fixed angle relative to the second, 
			and a fixed dihedral relative to the first two
	note "dihedral" really means torsion. this is an abstract 3D random walk
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

	#---destination
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

	#---! topology example (remove when finished)
	#topex = GMXTopology('inputs/test.top')

	#---note that the on-lattice version had a custom topology writer. here we try the automatic way.
	top = GMXTopology()

	#---construct the new molecule
	new_mol = {'moleculetype':{'molname':'AGLC','nrexcl':3},
		'atoms':[],'bonds':[],'angles':[],'dihedrals':[]}

	#---!
	mol = GMXTopology(state.aglc_source)
	#---select the correct molecule from the ITP
	molspec = mol.molecules['Other']

	#---!!! THE QTOT IS WRONG AND WHEN YOU FIX THE TOP FILE YOU SHOULD FIX IT!!!!!

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

	#---destination
	cwd = state.here
	n_p = state.melt_settings['n_p']

	#---instantiate a monomer to retrieve the correct spacing of the repeat unit (a0)
	#mono = Monomer(gro='aglc-cg.gro',remove_rules=False)

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
	history of manipulating which bonds get included:
		2017.06.15
			started with a previous run that only had sections 1-4
				and the backbone was pretty straight 
				but the sidechains were wobbly (see "arch-v024-before-adding-124-angles")
			while reviewing the method, we added section 5 
				and noticed that the sidechains were less wobbly
				for some reason the dihedrals of the sidechains looked really straight
				but there were no dihedrals in the itp file. weird.
			added section 6 and made the torsion for the backbone 142.0 from the melts_tuner.py
	"""

	#---!!!
	if False:

		#---! HANDLE THE LAST ONE MANUALLY
		mnum = n_p-1
		new_mol['angles'].append(dict(angle_template_flank_ortho,i=(mnum+1)*3+1+1,j=mnum*3+1+0,k=(mnum+1)*3+1))
		new_mol['angles'].append(dict(angle_template_flank_ortho,i=(mnum+1)*3+1+2,j=mnum*3+1+0,k=(mnum+1)*3+1))

	top.add_molecule(**{polymer_molame:new_mol})
	top.write(state.here+'dextran.itp')
	state.itp = ['dextran.itp']
	component(polymer_molame,count=1)

###---ADVANCED METHOD (started Nov 2017)

def planeproject(x,n): return x-np.dot(x,n)/np.linalg.norm(n)*vecnorm(n)

def vecangle(v1,v2):
	"""
	Compute the anglebetween two vectors
	"""
	v1n,v2n = vecnorm(v1),vecnorm(v2)
	dotp = np.dot(v1n,v2n)
	angle = np.arccos(dotp)*(180./np.pi)	
	if np.isnan(angle): return (0.0 if (v1n==v2n).all() else np.pi*(180/np.pi))
	return angle

def measure_torsions(coords,subsel):
	"""
	???
	"""
	dihedrals_measured = []
	# loop over instances of a particular bond
	for inds in subsel:
		# functions for computing torsions
		bulk_norm = lambda x: (x/np.tile(np.linalg.norm(x,axis=1),(3,1)).T)
		# project the AB vector onto the BC plane
		vecsBC = coords[:,inds[3]]-coords[:,inds[2]]
		vecsBCn = bulk_norm(vecsBC)
		vecsBAn = bulk_norm(coords[:,inds[2]]-coords[:,inds[1]])
		#! bulk plane project via:
		#! ... def planeproject(x,n): return x-np.dot(x,n)/np.linalg.norm(n)*vecnorm(n)
		vecsBA_projected_BC = (
			vecsBAn-np.tile(np.array([np.dot(m,n) 
				for m,n in zip(vecsBAn,vecsBCn)]),(3,1)).T*vecsBCn)
		# project the AB vector onto the BC plane
		vecsCB = coords[:,inds[1]]-coords[:,inds[2]]
		vecsCBn = bulk_norm(vecsCB)
		vecsCDn = bulk_norm(coords[:,inds[3]]-coords[:,inds[2]])
		vecsCD_projected_CB = (
			vecsCDn-np.tile(np.array([np.dot(m,n) 
				for m,n in zip(vecsCDn,vecsCBn)]),(3,1)).T*vecsCBn)
		# take the angle difference between the two projected vectors
		dihedrals_measured.append(np.array([vecangle(m,n) for m,n in zip(vecsBA_projected_BC,vecsCD_projected_CB)]))
	return np.array(dihedrals_measured)

class MultiScaleModel:
	"""
	Manage the mappings between scales.
	"""
	def __init__(self,**kwargs):
		"""Unpack a mapping into a model."""
		mapping = kwargs.pop('mapping')
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		style = mapping.pop('style',None)
		#---mappings which supply different "parts" (i.e. terminus or middle) of a polymer
		if style == 'coarse_to_atomistic': 
			self.style = style
			self.parts = mapping.pop('by_part',{})
		else: raise Exception('invalid style %s'%style)
		if mapping: raise Exception('unprocessed mapping directives %s'%mapping)
		#---shared variables
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
	def imagine(self):
		"""
		"""
		import ipdb;ipdb.set_trace()
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
		#---collect the intra-residue bonds
		bonds_intra,bond_parts = [],[]
		#---! this is currently hard-coded. check the style later
		model_by_part = model_spec['model']['by_part']
		#---! insert a check to make sure the model has all the parts
		#---! save partnames to the datastructure?
		partnames = list(set([i['part'] for i in rep['beads']]))

		#---make a blank GMXTopology (note that you must always include the name here)
		new_top = GMXTopology()

		#---per-atom topology column requirements to be read from the model spec
		requirements = {
			'atoms':['charge','mass','type'],}
		#---quicker lookups between residues
		def get_bead_in_residue(bead_name,resnum,beads_name='beads'):
			matches = [ii for ii,i in enumerate(rep[beads_name]) 
				if i['name']==bead_name and i['resnum']==resnum]
			if len(matches)>1: raise Exception
			elif len(matches)==0: return None
			else: return matches[0]
		#---construct the atoms list which is required for indexing the bonds later on
		atom_id,atoms = 1,[]
		bonds,angles,dihedrals,bonds_abstract = [],[],[],[]
		#---loop over residues
		for resnum in rep['resnums_u']:
			#---get the beads in this residue
			beads_inds,beads = zip(*[(ii,i) for ii,i in enumerate(rep['beads']) if i['resnum']==resnum])
			#---! first bead is the resname
			resname = beads[0]['part']
			#---loop over beads
			for atom_residue_index,(index,bead) in enumerate(zip(beads_inds,beads)):
				#---construct an atom entry
				#---several columns are automatic or identical to the abstract model
				new_atom = dict(id=atom_id,resnr=resnum,cgnr=resnum)
				atom_id += 1
				#---! is this elegant?
				new_atom.update(resname=model_spec['model']['by_part'][resname]['resname'])
				new_atom.update(atom=bead['name'])
				#---collect atom-specific columns from the model spec
				for key in requirements['atoms']:
					#---! ryan is worried this is still too complicated?
					new_atom[key] = model_spec['model']['by_part'][resname][
						'atoms'][atom_residue_index][bead['name']][key]
				atoms.append(new_atom)
			#! moving the following below to the automatic part
			if False:
				#---loop over bonds
				for bead_1,bead_2 in model_spec['model']['by_part'][resname].get('bonds',[]):
					lmatch = get_bead_in_residue(bead_name=bead_1,resnum=resnum)
					rmatch = get_bead_in_residue(bead_name=bead_2,resnum=resnum)
					bead_1_i,bead_2_i = lmatch,rmatch
					#---generate a new bond
					#---! inferring the indices here, starting at one instead of zero (use atom_id instead?)
					new_bond = dict(i=bead_1_i+1,j=bead_2_i+1)
					new_bond['funct'] = 1
					new_bond['force'] = 10000
					if not general_bonds:
						#---get the positions for this bond pair
						positions = coords[:,np.array([bead_1_i,bead_2_i])]
						#---compute observed distances
						distances = np.linalg.norm(positions.ptp(axis=1),axis=1)
						new_bond['length'] = distances.mean()
					else:
						# previously we got the exact index on coords (reduced i.e. coarse-grained positions) 
						# ... from the atomistic simulation. now we read the coords by name to generalize it
						# get all rows in coods that match the name of bead_1 and bead_2
						lmatch,rmatch = [np.array([i for i in [get_bead_in_residue(
							bead_name=bead,resnum=r,beads_name='beads_base') 
							for r in rep.get('resnums_u_base',rep['resnums_u'])] if i!=None]) 
							for bead in [bead_1,bead_2]]
						try: distances = (coords[:,lmatch]-coords[:,rmatch])
						except:
							import ipdb;ipdb.set_trace()
						#! note that distances is nframes by the number of redudant beads of that type by XYZ
						distances_mean = np.array([np.linalg.norm(bd.T,axis=1).mean() 
							for bd in distances.T]).mean()
						new_bond['length'] = distances_mean
					bonds.append(new_bond)
			# get connectivity here
			for bead_1,bead_2 in model_spec['model']['by_part'][resname].get('bonds',[]):
				lmatch = get_bead_in_residue(bead_name=bead_1,resnum=resnum)
				rmatch = get_bead_in_residue(bead_name=bead_2,resnum=resnum)
				bead_1_i,bead_2_i = lmatch,rmatch
				new_bond = dict(i=bead_1_i+1,j=bead_2_i+1)			
				bonds_abstract.append(new_bond)
		#---handle bonds between residues
		resnums_adjacent = [(rep['resnums_u'][i],rep['resnums_u'][i+1]) 
			for i in range(len(rep['resnums_u'])-1)]
		between_dict = model_spec['model'].get('between_parts',{})
		for left,lspec in between_dict.items():
			for right,rspec in lspec.items():
				#---search for all bonds between adjacent residues
				for bead_1,bead_2 in rspec.get('bonds',[]):
					#---check all adjacent residues
					for r1,r2 in resnums_adjacent:
						lmatch = get_bead_in_residue(bead_name=bead_1,resnum=r1)
						rmatch = get_bead_in_residue(bead_name=bead_2,resnum=r2)
						if lmatch and rmatch: 
							# previously we were calculating bonds here but this is moved below
							if False:
								#! this section is repetitive with the bond creation routine above
								bead_1_i,bead_2_i = lmatch,rmatch
								#---get the positions for this bond pair
								positions = coords[:,np.array([bead_1_i,bead_2_i])]
								#---compute observed distances
								distances = np.linalg.norm(positions.ptp(axis=1),axis=1)
								#---generate a new bond
								#---! inferring the indices here (use atom_id instead?)
								new_bond = dict(i=bead_1_i+1,j=bead_2_i+1)
								new_bond['funct'] = 1
								new_bond['force'] = 1250
								new_bond['length'] = distances.mean()
								bonds.append(new_bond)
							bead_1_i,bead_2_i = lmatch,rmatch
							new_bond = dict(i=bead_1_i+1,j=bead_2_i+1)
							bonds_abstract.append(new_bond)

		# automatically detect all possible angles 
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
		# loop over angles and dihedrals
		for nbeads in [2,3,4]:
			bead_candidates = list(set([tuple(sorted([i.num for i in a])) 
				for a in get_paths(molecule.atoms,d=nbeads)]))
			valid_bonds = [i for i in bead_candidates if len(set(list(i)))==nbeads]
			# loop over these bonds
			for beads in sorted(valid_bonds):
				if nbeads==2:
					i,j = beads
					inds = [m-1 for m in beads]
					new_bond = dict(i=i,j=j)
					if not general_bonds:
						#---get the positions for this bond pair
						positions = coords[:,np.array(inds)]
						#---compute observed distances
						distances = np.linalg.norm(positions.ptp(axis=1),axis=1)
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
						distances_mean = np.array([np.linalg.norm(bd.T,axis=1).mean() 
							for bd in distances.T]).mean()
						new_bond['length'] = distances_mean
					new_bond['funct'] = 1
					new_bond['force'] = 1000
					bonds.append(new_bond)
				elif nbeads==3:
					#! beware python 2 will reset i in a comprehension hence this is a poor index choice
					i,j,k = beads
					#! back to zero indexes
					inds = [m-1 for m in beads]
					new_angle = dict(i=i,j=j,k=k)
					new_angle['funct'] = 2
					new_angle['force'] = 500. # 2000 is too high for all of the bonds. need to select them
					# measure the angles 1-1 from the underlying coordinates
					if not general_bonds:
						vecs1,vecs2 = np.array((
							coords[:,inds[0]]-coords[:,inds[1]],
							coords[:,inds[2]]-coords[:,inds[0]]))
						angles_measured = np.array([vecangle(i,j) for i,j in zip(vecs1,vecs2)])
						new_angle['angle'] = angles_measured.mean()
					# the coarse model may have more coordinates than the (reduced) source
					else: 
						# to get the generic angles we need to check the inds list against the valid bond list
						# ... because just checking by name means we end up looking for all of the SB2, then
						# ... all of the SB3, then all of the MB1 and we get 1,1,3 candidates. instead of this
						# ... method we need to generalize the entire angle at once.
						bead_names = [rep['beads'][b]['name'] for b in inds]
						#! vital note: we prepared to use resnums to ensure that the bead names were also 
						#! ... split across residue numbers correctly. however, this was not needed since the #! ... following  section only selects triplets from the reference structure and 
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
							for i,j in zip(vecsAB[ii],vecsCB[ii])]) for ii in range(len(vecsAB))])
						# average over instances and frames
						new_angle['angle'] = angles_measured.mean()
					angles.append(new_angle)
				elif nbeads==4:
					if not general_bonds:
						# under development
						#! beware python 2 will reset i in a comprehension hence this is a poor index choice
						i,j,k,l = beads
						#! back to zero indexes
						inds = [m-1 for m in beads]
						new_dihedral = dict(i=i,j=j,k=k,l=l)
						new_dihedral['funct'] = 1
						new_dihedral['multiplicity'] = 1 # MAY NOT REALLY BE ONE!
						new_dihedral['force'] = 50.
						# COMPUTE TORSIONS! CHECK THIS LATER!
						dihedrals_measured = measure_torsions(coords,[inds])
						new_dihedral['angle'] = dihedrals_measured.mean()
						dihedrals.append(new_dihedral)
					else:
						# mimic the angle section
						inds = [m-1 for m in beads]
						i,j,k,l = beads 
						new_dihedral = dict(i=i,j=j,k=k,l=l)
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
						new_dihedral['force'] = 50.
						new_dihedral['angle'] = dihedrals_measured.mean()
						dihedrals.append(new_dihedral)
				else: raise Exception('development note. cannot generate abond with %d beads'%nbeads)

		if remember_bonds:
			self.reps[source]['bond_info'] = dict(bonds=bonds,angles=angles,dihedrals=dihedrals)
		else:
			#---populate the topology
			#---! note that you have to pass in moleculetype or the ITP is incomplete
			new_top.molecules['DEX'] = dict(moleculetype=dict(molname='DEX',nrexcl=3),
				atoms=atoms,bonds=bonds,angles=angles,dihedrals=dihedrals)
			new_top.write(state.here+'dextran.itp')

		"""
		note that in the above we calculate distances, angles, etc using the model positions
		which assumes the CG model has the same number of residues as the atomistic model
		we should obviously relax this requirement
		"""

def make_polymer(name='melt',**kwargs):
	"""
	Make a melt by placing monomers according to an off-lattice walk.
	This is the first section of make_cg_gel; the second section built the topology.
	"""
	polymer_molname = 'DEX'
	monomer_resname = 'AGLC'

	#---destination
	cwd = state.here
	n_p = state.melt_settings['n_p']
	angle,torsion = state.melt_settings['angle'],state.melt_settings['torsion']
	a0 = state.melt_settings['a0']
	#---prepare an abstract 3D walk
	walk_abstract_pts = make_off_lattice_walk(n_p-1,a0,angle,torsion)
	points_with_sides = []
	for pt in walk_abstract_pts:
		points_with_sides.append(pt)
		for i in range(2):
			#---! randomly place sidechain points 1 a0 distance away. map 0 to 1 to -1 to 1
			points_with_sides.append(pt+(1-2*np.random.rand(3))/10.)
	points_with_sides = np.array(points_with_sides)
	residue_names = np.array([monomer_resname for p in points_with_sides])
	residue_indices = np.array([i/3+1 for i in range((n_p)*3)])
	#! hardcoding the naming scheme here but it would be useful to get this from the YAML file!
	atom_names = np.array(['%sB%d'%(l,i) for l in 'S'+'M'*(n_p-2)+'E' for i in range(1,3+1)])
	polymer = GMXStructure(pts=points_with_sides,residue_indices=residue_indices,
		residue_names=residue_names,atom_names=atom_names,box=[10.,10.,10.])
	polymer.write(state.here+'%s-built.gro'%name)
	gmx('editconf',structure='%s-built'%name,gro=name,c=True,log='editconf-center-polymer')
	#---add to the composition
	state.itp = ['dextran.itp']
	component(polymer_molname,count=1)

def make_crude_coarse_polymer(name='vacuum',diagonals=False,review=False):
	"""
	Settings in the experiment bring you here, or above to `make_polymer`.
	The `make_polymer` code above will make a single large molecule. 
	This method will make a melt from a square lattice. 
	UNDER DEVELOPMENT.
	"""
	monomer_resname = 'AGLC'
	polymer_molname = 'DEX'
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
		while size<n_p-1:
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
			if pnum==0: atom_name_letter = 'S'
			elif pnum<n_p-1: atom_name_letter = 'M'
			elif pnum==n_p-1: atom_name_letter = 'E'
			else: raise Exception
			atom_names.extend(['%sB%d'%(atom_name_letter,i) for i in range(1,3+1)])
			residue_indices.extend([wnum+1 for i in range(3)])
	coords = np.array(coords)*a0
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
	state.itp = ['dextran.itp']

def forward_mapper():
	"""
	Compute statistics from an atomistic polymer on a coarse-grained mapping.
	Used to perform a "forward mapping" from atomistic to coarse-grained.
	"""
	#---get the reference structure, presumably from another run
	ref_specs = settings.atomistic_reference
	ref_gro = os.path.join(ref_specs['path'],ref_specs['gro'])
	ref_xtc = os.path.join(ref_specs['path'],ref_specs['xtc'])
	uni = MDAnalysis.Universe(ref_gro,ref_xtc)
	sel = uni.select_atoms(ref_specs['selection'])
	
	#---read the mapping
	import yaml
	with open(settings.mapping_spec) as fp: mapping = yaml.load(fp.read())
	model = MultiScaleModel(**mapping)

	#---! move this into the class?
	lenscale = 10.0
	mass_table = {'H':1.008,'C':12.011,'O':15.999,'N':14.007,'P':30.974}
	#---the atomistic reference tells us how to get molecules from that simulation
	molecule_map = ref_specs.get('molecule_map',None)
	#---read the statistics for a single molecule
	#---prepare a crude structure for the atomistic points
	#---! this block might be included in the model class but for now it's only tangent
	fine = dict(molecule_map=molecule_map,selection=ref_specs['selection'])
	fine.update(resnums=sel.resnums,names=sel.names)
	fine.update(masses=np.array([mass_table[i[0]] for i in fine['names']]))
	model.add_representation('fine',**fine)
	# after preparing fine we make the coarse model
	if molecule_map == 'one molecule many residues injective':
		model.build_fine_to_coarse('fine','coarse')
	elif molecule_map == 'one molecule many residues bigger':
		model.build_fine_to_coarse('fine','coarse',molecule_map='one molecule many residues injective')
		specs = dict(nres=settings.melt_settings['n_p'],molecule_map='one molecule many residues bigger')
		model.build_fine_to_coarse('fine','coarse_big',**specs)
	else: raise Exception('invalid molecule map %s'%molecule_map)

	#---reduce each frame of theectory to coarse-grained coordinates
	nframes = len(uni.trajectory)
	frame_skip = 100
	n_cg_beads = model.reps['coarse']['length']
	coords_red = np.zeros((len(range(0,nframes,frame_skip)),n_cg_beads,3))
	for ff,fr in enumerate(range(0,nframes,frame_skip)):
		uni.trajectory[fr]
		#---apply Angstroms to nm conversion here
		coords_red[ff] = model.reduce(source='fine',target='coarse',positions=sel.positions/lenscale)

	#---! crude method to write the coarse-grained trajectory for review
	if True:
		residue_indices = model.reps['coarse']['resnums']
		atom_names = model.reps['coarse']['names']
		residue_names = np.array(['AGC' for i in residue_indices])
		for fr in range(len(coords_red)):
			polymer = GMXStructure(pts=coords_red[fr],residue_indices=residue_indices,
				residue_names=residue_names,atom_names=atom_names,box=[100.,100.,100.])
			polymer.write(state.here+'cg_model%04d.gro'%fr)
		bash('cat cg_model0* > cg_model.gro',cwd=state.here)
		"""
		You can review the mapping by looking at the trajectory dumped above.
		In the original version of this code, we took another frame and tried backmapping.
		See the experiments file for more development notes.
		"""

	#---! inspecting the distributions (note that this method was folded into interpret_model)
	if False:
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

	#---read the mapping
	with open(settings.model_spec) as fp: model_spec = yaml.load(fp.read())
	#---create a model from the mapping
	#! make the following line contingent on whether we are doing injective or bigger?
	model.interpret_model('coarse',coords_red,**model_spec)
	model.interpret_model('coarse_big',coords_red,**model_spec)
