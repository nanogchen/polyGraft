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

import MDAnalysis as mda
import numpy as np

class Polymer():

	def __init__(self, poly_name):
		"""atoms, n_atoms, bonds, angles, dihedrals
		atoms.names/types/resnames/resids/segids/positions

		"""
		self.polyname_ = poly_name
		self.polyGRO_ = None
		self.polyITP_ = None		

	def readGRO(self, GROfile):
		# read from files
		self.polyGRO_ = mda.Universe(GROfile)

		# center the polymer
		self.center2graft()

	def readITP(self, ITPfile):
		# read from files
		self.polyITP_ = mda.Universe(ITPfile, topology_format="ITP")

	def getOrient(self):

		""" get the principal axis of a polymer chain for grafting use; 
		when polymer is aligned in a direction, 
		then the component in that direction is the smallest. 
		then the ref direction is the third principal axis, 
		the arrow is point from sulfur to the com"""
		
		_,_,e3 = self.polyGRO_.atoms.principal_axes()

		r1 = self.polyGRO_.atoms.positions[0]
		com = self.polyGRO_.atoms.center_of_geometry()
		vect = com - r1

		if np.dot(e3, vect) > 0:
			return e3
		else:
			return e3*(-1)

	def getEndVect(self):
		"""
		get the end-to-end vector of the chain from first and last index
		"""
		ref_vec = self.polyGRO_.atoms.positions[-1] - self.polyGRO_.atoms.positions[0]
		ref_vec = ref_vec/np.linalg.norm(ref_vec)
		return ref_vec

	def center2graft(self):
		"""
		move the grafting point (first atom) of a polymer to the origin 
		"""

		pos = self.polyGRO_.atoms.positions
		self.polyGRO_.atoms.positions = pos - self.polyGRO_.atoms.positions[0]

	def align(self, principal_axes_idx, ref_direction):
		"""align one principal axis of a polymer to a given direction """

		assert principal_axes_idx in [0,1,2], f"Only 0/1/2 of the axis is allowed for align the polymer!"

		self.polyGRO_.atoms.align_principal_axis(principal_axes_idx, ref_direction)
		