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
from numpy.linalg import norm
from numba import njit, prange
import math

@njit
def getNN(xyz,nearest_neighbor):
	"""
	get the nearest bond info of given atoms
	"""

	# nearest dist
	tol = nearest_neighbor*0.1

	# bonds
	Natoms = xyz.shape[0]
	Max = Natoms*(Natoms-1)
	bonds = np.empty((Max, 2), dtype=np.int_)	

	# loop
	nbonds = 0
	for i in prange(Natoms):
		for j in prange(i+1,Natoms):

			vec_ij = xyz[i] - xyz[j]

			if norm(vec_ij) < nearest_neighbor + tol:
				bonds[nbonds,:] = np.asarray([i,j])
				nbonds += 1						

	return bonds

class Crystal():

	def __init__(self, *args):
		assert len(args) > 1, f"Constructor of Crystal should specify: type, element, and parameter"
		assert args[0] in ['nanorod', 'nanoparticle','nanopore','nanoslab'], f"lattice_type can ONLY be one of nanorod/nanoparticle/nanopore/nanoslab, but {lattice_type} was given!"

		if args[0] == "nanoparticle":
			assert len(args) == 3, f"For nanoparticle (1), only element (2) and radius (3) is needed!"
			self.pars_ = [args[2]]

		elif args[0] == "nanopore" or args[0] == "nanorod":
			assert len(args) == 4, f"For nanopore/nanorod (1), only element (2) radius (3) and depth (4) are needed!"
			self.pars_ = [args[2], args[3]]

		elif args[0] == "nanoslab":
			assert len(args) == 5, f"For nanoslab (1), only element (2) length (3), width (4) and depth (5) are needed!"
			self.pars_ = [args[2], args[3], args[4]]

		# this will be a Universe object
		self.crystal_ = None
		self.crystal_type_ = args[0]
		self.crystal_element_ = args[1]
	
	def readDATA(self, fname, atom_style="id type x y z", guessing_bond=False, nearest_neighbor=1):
		# read lammps data file

		# set crystal_ based on atomsk generated structure
		self.crystal_ = mda.Universe(fname, atom_style=atom_style)
		self.setAtomTypes()
		self.setMass()
		self.setCharge(0.0)

		# generate nearest neighbors
		if guessing_bond:

			# bonds
			xyz = self.crystal_.atoms.positions
			bonds = getNN(xyz,nearest_neighbor)

			# find effective bonding
			bonds_nn = bonds[~np.all(bonds == 0, axis=1)]

			self.crystal_.add_TopologyAttr('bonds', bonds_nn)
			
	def setAtomTypes(self):
		atom_types = [self.crystal_element_ for _ in range(self.crystal_.atoms.n_atoms)]

		self.crystal_.add_TopologyAttr('type', atom_types)
		self.crystal_.add_TopologyAttr('names', atom_types)
		self.crystal_.add_TopologyAttr('resnames', ["LAT"])

	def setMass(self):
		"""set the atom mass of the crystal"""
		atom_types = [self.crystal_element_ for _ in range(self.crystal_.atoms.n_atoms)]
		mass = mda.topology.guessers.guess_masses(atom_types)	
		self.crystal_.add_TopologyAttr('mass', mass)

	def setCharge(self, charge=0.0):
		"""set the atom mass of the crystal"""
		atom_charges = [charge for _ in range(self.crystal_.atoms.n_atoms)]
		self.crystal_.add_TopologyAttr('charge', atom_charges)




