#!/usr/bin/env python

_not_reported = ['vecnorm','dotplace','rotation_matrix','check_cols','review3d','apply_rotation','align']
#---! somehow status was getting reported when I merged in melts_tuner.py
#_not_reported += ['status']

#vecnorm = lambda x: x/np.linalg.norm(x)
#dotplace = lambda n: re.compile(r'(\d)0+$').sub(r'\1',"%8.3f"%float(n)).ljust(8)

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
