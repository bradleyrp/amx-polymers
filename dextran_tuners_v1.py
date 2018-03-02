#!/usr/bin/env python

"""
convert normal distribution standard deviations to force constants. note that this
	also serves as a filter to exclude certain bonds with wide distributions. see the 
	bonds_review_width figures to see where to draw the cutoff
"""

def filter_angles(std,**kwargs):
	#! previously we stratified angles by the gaussian widths
	return 0.
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
