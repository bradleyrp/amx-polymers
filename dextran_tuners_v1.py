#!/usr/bin/env python

"""
convert normal distribution standard deviations to force constants. note that this
	also serves as a filter to exclude certain bonds with wide distributions. see the 
	bonds_review_width figures to see where to draw the cutoff
"""

def filter_angles_v1(std,**kwargs):
	#! previously we stratified angles by the gaussian widths
	return 0.
	#! example
	if False:
		if std<=3: return 25.
		else: return 0.
	#! strong angles inside monomers
	if set(kwargs['names']) in [{'SB1','SB2','SB3'},{'MB1','MB2','MB3'},{'EB1','EB2','EB3'}]:
		return 100.
	else: return 20.

def filter_dihedrals(std,**kwargs):
	return 0 #! dihedrals cause crashes!
	if std<=4: return 5.
	else: return 0.

def filter_bonds(std,**kwargs): 
	return 10000.

def filter_angles(std,**kwargs):
	#! be very careful to retain symmetry here
	bond_table = {
		('MB3','MB1','MB2'):20,('MB2','MB1','MB3'):20,
		('MB1','MB2','MB3'):50,('MB3','MB2','MB1'):50,}
	return bond_table.get(tuple(kwargs['names']),0)
