#!/usr/bin/env python

"""
convert normal distribution standard deviations to force constants. note that this
	also serves as a filter to exclude certain bonds with wide distributions. see the 
	bonds_review_width figures to see where to draw the cutoff
"""

def filter_dihedrals(std,**kwargs):
	return 0
	bond_table = {
		('MB1','MB2','MB3','MB2'):5,('MB2','MB3','MB2','MB1'):5,}
	return bond_table.get(tuple(kwargs['names']),0)

def filter_bonds(std,**kwargs): 
	return 30000.

def filter_angles(std,**kwargs):
	names = tuple(kwargs.pop('names'))
	resids = kwargs.pop('resids',None)
	if not resids: raise Exception
	bond_mappings = [
		{'names':('MB1','MB2','MB3'),'resids':'000','strength':25},
		{'names':('MB1','MB2','MB1'),'resids':'001','strength':25},
		{'names':('MB3','MB2','MB1'),'resids':'001','strength':25},
		{'names':('MB2','MB1','MB2'),'resids':'011','strength':25},]

	matches = []
	for bond in bond_mappings:
		if bond['names']==names or bond['names'][::-1]==names:
			if bond['resids'] == ''.join([str(i-min(resids)) for i in sorted(resids)]):
				matches.append(bond['strength'])
	#import ipdb;ipdb.set_trace()
	if len(matches)==0: return 0
	elif len(matches)==1: return matches[0]
	else: raise Exception('redundant bond list matches: %s'%matches)

def scale_bonds(l):
	bot,top,inc = 0.26,0.56,0.02
	return min([round((max([l-bot,0]))/inc)*inc+bot,top])
