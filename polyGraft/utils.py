# 
# Copyright (C) Guang Chen
# 
# This file is part of polyGraft program
#
# polyGraft is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# polyGraft is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
""" custom functions called from the other classes"""

import numpy as np
import math
import sys

def getTransformationMat(ref_vec, norm_vector):
	""" get the transformation matrix using the normal vector at a grafting point 
	map from the local ref_vect to the global norm_vector
	"""

	# Au-S bond length 2.65 A
	Au_S = 2.65

	# ref set in x direction
	# ref_vec = np.asarray([1.0, 0.0, 0.0])
	ref_vec = ref_vec/np.linalg.norm(ref_vec)
	norm_vector = norm_vector/np.linalg.norm(norm_vector)

	vect_k = np.cross(ref_vec, norm_vector)	
	vect_u = vect_k/np.linalg.norm(vect_k) # have to normalize the rotation axis (unit vector)
	skew_mat = get_skew_matrix(vect_u)
	sin_theta = np.linalg.norm(vect_k)/(np.linalg.norm(ref_vec)*np.linalg.norm(norm_vector))
	cos_theta = np.dot(ref_vec, norm_vector)/(np.linalg.norm(ref_vec)*np.linalg.norm(norm_vector))

	RotMat = np.identity(3) + skew_mat*sin_theta + np.matmul(skew_mat,skew_mat)*(1-cos_theta) 

	# assert a rotation matrix
	assert (np.linalg.det(RotMat) - 1) < 1e-5, f"The rotation matrix should have determinant of 1, but got {np.linalg.det(RotMat)}"

	# translational matrix
	TransMat = norm_vector*Au_S
	TransMat = np.expand_dims(TransMat, axis=0)

	return RotMat,TransMat

def get_skew_matrix(vector):
	"""get the skew matrix of a given vector
	"""

	mat = np.zeros((3,3), dtype=float)
	mat[0,1] = -vector[2]
	mat[0,2] = vector[1]
	mat[1,0] = vector[2]
	mat[1,2] = -vector[0]
	mat[2,0] = -vector[1]
	mat[2,1] = vector[0]

	return mat

def getNN_two(pos1, pos2):

	"""
	get the nearest neighbor (atom2) of atom 1
	"""

	atom_idx = []

	# loop each sulfur atom
	for ipos in pos1:

		# np.argmin
		dist_vec = np.linalg.norm(pos2-ipos, axis=1)
		gft_idx = np.argmin(dist_vec)
		atom_idx.append(gft_idx)

	return pos2[atom_idx,:], np.array(atom_idx)+1

def rad2deg(rad):

	return 180.0*rad/math.pi

def deg2rad(deg):

	return deg*math.pi/180.0

def getUniqueBondTypes(mol):

	bondtype = []
	for ibond in mol.atoms.bonds:
		bondtype.append(ibond.type)

	return list(set(bondtype))	

def getUniqueAngleTypes(mol):

	angletype = []
	for iangle in mol.atoms.angles:
		angletype.append(iangle.type)

	return list(set(angletype))

def getUniqueDihedralTypes(mol):

	dihedraltype = []
	for idihedral in mol.atoms.dihedrals:
		dihedraltype.append(idihedral.type)

	return list(set(dihedraltype))

def getUniqueImproperTypes(mol):

	impropertype = []
	for iimproper in mol.atoms.impropers:
		impropertype.append(iimproper.type)

	return list(set(impropertype))