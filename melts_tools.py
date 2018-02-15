#!/usr/bin/env python

import numpy as np
import re

### GEOMETRY

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

### VISUALIZATION

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

# additions for use in MultiScaleModel

def planeproject(x,n): 
	"""Project points `x` onto a plane given by a normal vector `n`."""
	return x-np.dot(x,n)/np.linalg.norm(n)*vecnorm(n)

def vecangle(v1,v2):
	"""Compute the anglebetween two vectors"""
	v1n,v2n = vecnorm(v1),vecnorm(v2)
	dotp = np.dot(v1n,v2n)
	angle = np.arccos(dotp)*(180./np.pi)	
	if np.isnan(angle): return (0.0 if (v1n==v2n).all() else np.pi*(180/np.pi))
	return angle

def measure_torsions(coords,subsel):
	"""Measure torsions from a sub-selection of some coordinates."""
	dihedrals_measured = []
	# loop over instances of a particular bond
	for inds in subsel:
		# functions for computing torsions
		bulk_norm = lambda x: (x/np.tile(np.linalg.norm(x,axis=1),(3,1)).T)
		# project the AB vector onto the BC plane
		vecsBC = coords[:,inds[3]]-coords[:,inds[2]]
		vecsBCn = bulk_norm(vecsBC)
		vecsBAn = bulk_norm(coords[:,inds[2]]-coords[:,inds[1]])
		# bulk plane project via:
		# ... def planeproject(x,n): return x-np.dot(x,n)/np.linalg.norm(n)*vecnorm(n)
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
