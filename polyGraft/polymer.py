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

# bond/angle parameters from opls-aa
pars_dict = {
	"CC" : 1.529,
	"CO" : 1.41,
	"OC" : 1.41,
	"OCC" : 109.5,
	"COC" : 109.5,
	"CCO" : 109.5,
	"CCC" : 112.7
}

from utils import rad2deg,deg2rad
import MDAnalysis as mda
import numpy as np
import math
import sys

def getRadius(Npva):
	"""get the radius of the cyclic backbone"""

	ccbond = pars_dict['CC']
	rad = 0.5*ccbond/math.sin(0.5*math.pi/Npva)

	return rad

class Polymer():

	def __init__(self, poly_name):
		"""atoms, n_atoms, bonds, angles, dihedrals
		atoms.names/types/resnames/resids/segids/positions

		"""
		self.polyname_ = poly_name
		self.Nrepeats_ = 0
		self.topology_ = None
		self.polyGRO_ = None
		self.polyITP_ = None

	def setNrepeats(self, Nrepeats):
		self.Nrepeats_ = Nrepeats

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

	def updatePos(self, pos):
		assert pos.shape[0] == self.polyGRO_.atoms.n_atoms, f"The given position should have the same length as the NO. of atoms!"
		self.polyGRO_.atoms.positions = pos

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

	def gen_pos(self, Nrepeats, topology='linear'):
		# generate the position of the given polymer: heavy atoms only (hydrogen will be taken care of by Avogadro)

		pos = []
		atomnames = []
		self.Nrepeats_ = Nrepeats
		self.topology_ = topology

		if self.polyname_ == "PEO":

			if topology == 'linear':

				# the basic length
				dx = pars_dict['CC']*math.cos(deg2rad(90.0-0.5*pars_dict['CCC']))
				dy = pars_dict['CC']*math.sin(deg2rad(90.0-0.5*pars_dict['CCC']))

				for ires in range(Nrepeats):
					x_inc = ires*3*dx 

					if ires % 2 == 0:
						# first carbon
						xc = x_inc
						yc = 0.0
						zc = 0.0
						pos.append([xc,yc,zc])
						atomnames.append('C')

						# second carbon
						xc = x_inc+dx
						yc = -dy
						zc = 0.0
						pos.append([xc,yc,zc])
						atomnames.append('C')

						# oxygen
						xc = x_inc+2*dx
						yc = 0.0
						zc = 0.0
						pos.append([xc,yc,zc])
						atomnames.append('O')

					else:
						# first carbon
						xc = x_inc
						yc = -dy
						zc = 0.0
						pos.append([xc,yc,zc])
						atomnames.append('C')

						# second carbon
						xc = x_inc+dx
						yc = 0.0
						zc = 0.0
						pos.append([xc,yc,zc])
						atomnames.append('C')

						# oxygen
						xc = x_inc+2*dx
						yc = -dy
						zc = 0.0
						pos.append([xc,yc,zc])
						atomnames.append('O')
			else:
				print(f"Cyclic PEO is not available! Exiting...")
				sys.exit(0)

		elif self.polyname_ == "PVA":

			if topology == 'linear':
				# add head group
				pos.append([0.0, 0.0, 0.0])
				atomnames.append('C')

				# the basic length
				dx = pars_dict['CC']*math.cos(deg2rad(90.0-0.5*pars_dict['CCC']))
				dy = pars_dict['CC']*math.sin(deg2rad(90.0-0.5*pars_dict['CCC']))

				for ires in range(Nrepeats):
					x_inc = ires*2*dx

					# first carbon		
					xc = x_inc+dx
					yc = -dy
					zc = 0.0
					pos.append([xc,yc,zc])
					atomnames.append('C')

					# second carbon
					xc = x_inc+2*dx
					yc = 0.0
					zc = 0.0
					pos.append([xc,yc,zc])
					atomnames.append('C')

					# oxygen
					pos.append([xc,yc+pars_dict["OC"],zc])
					atomnames.append('O')

				# add tail
				pos.append([xc+dx, -dy, 0.0])
				atomnames.append('C')

			elif topology == 'cyclic':
				rad = getRadius(Nrepeats)
				theta = 2*math.pi/(2*Nrepeats)

				for iatom in range(Nrepeats):

					# first carbon		
					xc = rad*math.cos((2*iatom)*theta)
					yc = rad*math.sin((2*iatom)*theta)
					zc = 0.0
					pos.append([xc,yc,zc])
					atomnames.append('C')

					# second carbon
					xc = rad*math.cos((2*iatom+1)*theta)
					yc = rad*math.sin((2*iatom+1)*theta)
					zc = 0.0
					pos.append([xc,yc,zc])
					atomnames.append('C')

					# oxygen
					pos.append([xc,yc,zc+pars_dict["OC"]])
					atomnames.append('O')

			else:
				print(f"Unknown topology type ({topology}) other than linear and cyclic! Exiting...")
				sys.exit(0)

		else:
			print(f"Unknown polymer type ({self.polyname_}) other than PEO and PVA! Exiting...")
			sys.exit(0)

		# create an atom group
		self.polyGRO_ = mda.Universe.empty(n_atoms = len(pos))
		self.polyGRO_.add_TopologyAttr('names', atomnames)
		self.polyGRO_.add_TopologyAttr('elements', atomnames)
		self.polyGRO_.add_TopologyAttr('resnames', [self.polyname_])
		self.polyGRO_.load_new(np.array(pos), order='fac')
