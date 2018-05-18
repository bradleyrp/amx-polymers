#!/usr/bin/env python

#!!! the following functions taken from melts.py and adapted slightly. these need to be centralized
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

def get_paths(atoms,path=None,d=3):
	if d>0:
		if not path: path,next_up = [],atoms
		else: next_up = path[-1].neighbors
		for neighbor in next_up:
			local_path = path[:]+[neighbor]
			for p in get_paths(atoms,path=local_path,d=d-1): yield p
	else: yield path

def make_molecule(bonds_abstract):
	molecule = Molecule()
	bond_list = [tuple([bond[k] for k in 'ij']) for bond in bonds_abstract]
	for atom_id in list(set([i for j in bond_list for i in j])):
		molecule.atoms.append(Atom(atom_id))
	for bond in bond_list: molecule.add_link(bond)
	return molecule

import os

mapping = """
- atoms: [O2,C1,H12,H11,H13]
- atoms: [C3,C4,H4]
- atoms: [C5,C6,H5,H6]
- atoms: [C7,C8]
- atoms: [C9,O10,C11]
- atoms: [C42,C43,O44]
- atoms: [H131,O13,C12,C14]
- atoms: [H41,O41,C40,C15]
- atoms: [C20,H20,C19,H192,H193]
- atoms: [C16,C17,H162,H163,O18,H18]
- atoms: [O37,C36]
- atoms: [C38,H383,H382,O39,H39]
- atoms: [O22,C23,H23]
- atoms: [O35,C32,H32,C34,H341,H342,H343]
- atoms: [C29,H29,O31,H31]
- atoms: [C26,H26,N28,H281,H282,H283,C25,H252,H253]
"""

# ordered atom list
atoms_raw = "C3 C4 C5 C6 C7 C8 C11 C12 C14 C15 C40 C42 C9 O10 C16 C17 C19 C20 C43 O44 O22 C23 O35 C32 C29 C34 C26 C25 O31 N28 O13 O18 C36 O37 C38 O39 O41 O2 C1 H4 H5 H6 H162 H163 H192 H193 H20 H23 H32 H29 H26 H252 H253 H382 H383 H18 H39 H341 H342 H343 H31 H11 H12 H13 H281 H282 H283 H41 H131"

import yaml

bonds_abstract_yaml = """
- [1,2]
- [2,4]
- [2,6]
- [2,3]
- [3,4]
- [4,5]
- [4,6]
- [5,6]
- [5,7]
- [6,8]
- [5,8]
- [7,8]
- [7,10]
- [8,9]
- [7,9]
- [9,10]
- [10,11]
- [11,12]
- [9,13]
- [13,14]
- [13,16]
- [14,15]
- [15,16]
"""

def simple_coarse_grained_tutorial():
	"""
	"""
	maps = yaml.load(mapping)
	missing = set.difference(set(atoms_raw.split()),set([i for j in [k['atoms'] for k in maps] for i in j]))
	if missing: raise Exception('missing %s in the mapping'%missing)
	out = []
	for inum,item in enumerate(maps):
		out.append('[%d]'%(inum+1))
		out.append(' '.join(['%d'%(atoms_raw.split().index(i)+1) for i in item['atoms']]))
	with open(state.here+'mapping.ndx','w') as fp: fp.write('\n'.join(out)+'\n')
	n_beads_map = 16
	os.system('cd %s'%state.here+' && '+('seq 0 15 | gmx traj -f %(source_xtc)s -s %(source_gro)s -oxt '
		'mapped.gro -n mapping.ndx -com -ng %(n_beads_map)s')%dict(
		source_xtc=settings.source_xtc,source_gro=settings.source_gro,n_beads_map=n_beads_map))

	bonds_abstract = [{'i':i,'j':j} for i,j in yaml.load(bonds_abstract_yaml)]
	molecule = make_molecule(bonds_abstract)
	tops = []
	# loop over angles and dihedrals
	for nbeads in [2,3,4]:
		bead_candidates = []
		for a in get_paths(molecule.atoms,d=nbeads):
			this = [i.num for i in a]
			if (this not in bead_candidates and this[::-1] not in bead_candidates 
				and len(set(this))==nbeads):
				bead_candidates.append(this)
		valid_bonds = bead_candidates
		print(valid_bonds)
		tops.append('[%s]\n%s'%({2:'bonds',3:'angles',4:'dihedrals'}[nbeads],
			'\n'.join([' '.join([str(j) for j in i]) for i in valid_bonds]))+'\n')
	with open(state.here+'bonded.ndx','w') as fp: fp.write('\n'.join(tops))

	import ipdb;ipdb.set_trace()