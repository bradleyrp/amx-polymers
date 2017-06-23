#!/usr/bin/env python

import numpy as np
import MDAnalysis
import matplotlib as mpl
import matplotlib.pylab as plt

def get_angle_torsion_distn_from_pentamer():

	"""
	"""

	uni = MDAnalysis.Universe('../melt-v011/s01-melt/system.gro','../melt-v011/s01-melt/md.parts.pbcmol.xtc')
	sel = uni.select_atoms('resname AGLC')

	#---note that the string additions were added after Samaneh and Ryan sketched the key atom mapping
	#---...and the break below helped to show which atoms are 
	mapping = {
		'terminus_start':{
			'SB1':'C6 C5'+' H5 H61 H62 O6',
			'SB2':'O5 C4 C1 O1'+' H1 HO1 H4 O4 HO4',
			'SB3':'C3 C2'+' H2 O2 HO2 H3 O3 HO3'},
		'middle':{
			'MB1':'C6 C5'+' H5 H61 H62',
			'MB2':'O5 C4 C1 O6'+' H1 H4 O4 HO4',
			'MB3':'C3 C2'+' H2 O2 HO2 H3 O3 HO3'},
		'terminus_end':{
			'EB1':'C6 C5'+' H61 H62 HO6',
			'EB2':'O5 C4 C1 O6'+' H1 H5 H4 O4 HO4',
			'EB3':'C3 C2'+' H2 O2 HO2 H3 O3 HO3'}}

	#---! assume only a single AGLC molecule
	resnums = np.sort(np.unique(sel.atoms.resnums))
	term_start,term_end = resnums[0],resnums[-1]
	connector = lambda resnum: {0:'ERROR!',1:'terminus_start',len(resnums):'terminus_end'}.get(resnum,'middle') 

	#---iterate over residues to make sure we have mapped everything
	for resnum in np.unique(sel.atoms.resnums):
		#---figure out which connectivity we have 
		connectivity = connector(resnum)
		#---mapping happens by assigning indices to all of the atoms
		this_names = list(np.array(sel.residues[list(sel.residues.resnums).index(resnum)].atoms.names))
		this_map = dict([(k,np.array(v.split())) for k,v in mapping[connectivity].items()])
		for cg_bead_name,atoms in this_map.items():
			for atom in atoms: 
				if atom in this_names: this_names.remove(atom)
				else: print('%s not in %s'%(atom,this_names))
		#---report leftover names so you can add them to the mapping above
		if this_names: 
			print('breaking on %s'%connectivity)
			print(this_names)
			break

	#---map the resnum and atom name to the atom index
	map_resnum_name_index = dict([(i,ii) for ii,i in enumerate(zip(*[sel.atoms.resnums,sel.atoms.names]))])
	#---map residue number to connectivity
	map_resnum_connectivity = [connector(resnum) for resnum in sel.atoms.resnums]

	#---construct a list of resnums, bead_names, and corresponding atom indices
	cg_model = []
	for resnum in resnums:
		connectivity = connector(resnum)
		for bead_name in sorted(mapping[connectivity].keys()):
			atoms = mapping[connectivity][bead_name]
			atoms_in_bead = np.sort([map_resnum_name_index[(resnum,atom)] for atom in atoms.split()])
			cg_model.append([resnum,bead_name,atoms_in_bead])

	mass_table = {'H':1.008,'C':12.011,'O':15.999,'N':14.007,'P':30.974}

	nframes = len(uni.trajectory)
	frame_skip = 100
	coords_red = np.zeros((len(range(0,nframes,frame_skip)),len(cg_model),3))
	for ff,fr in enumerate(range(0,nframes,frame_skip)):
		for anum,(rnum,bead,atom_inds) in enumerate(cg_model):
			status('computing reduced coordinates',i=fr,looplen=nframes,tag='compute')
			uni.trajectory[fr]
			pos = sel.positions[atom_inds]/10.
			weights = np.array([mass_table[i[0]] for i in sel.names[atom_inds]])
			com = (pos*np.tile(weights,(3,1)).T/np.sum(weights)).mean(axis=0)
			#---! COG not COM below
			coords_red[ff][anum] = pos.mean(axis=0)

	residue_indices = list(zip(*cg_model)[0])
	atom_names = list(zip(*cg_model)[1])
	residue_names = np.array(['AGC' for i in residue_indices])
	for fr in range(len(coords_red)):
		polymer = GMXStructure(pts=coords_red[fr],residue_indices=residue_indices,
			residue_names=residue_names,atom_names=atom_names,box=[10.,10.,10.])
		polymer.write(state.here+'cg_model%04d.gro'%fr)

	bash('cat cg_model0* > cg_model.gro',cwd=state.here)

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

	#---measure the torsions on the backbone beads
	torsion_inds = np.array([[1,4,7,10],[4,7,10,13]])
	#---allocate angles
	torsion_angles = np.zeros((len(torsion_inds),len(coords_red)))
	for fr in range(len(coords_red)):
		status('computing torsions',i=fr,looplen=len(coords_red),tag='compute')
		for tt,torsion_ind in enumerate(torsion_inds):
			pts = coords_red[fr,torsion_ind]
			#---imagine the torsion is given by pts a,b,c,d
			#---project ba onto plane given by bc
			project_ba_bc = planeproject(pts[0]-pts[1],pts[2]-pts[1])
			project_cd_bc = planeproject(pts[3]-pts[2],pts[2]-pts[1])
			angle = vecangle(project_ba_bc,project_cd_bc)
			torsion_angles[tt,fr] = angle

	import ipdb;ipdb.set_trace()

	print('TORSION IS %.1f'%torsion_angles.mean())
