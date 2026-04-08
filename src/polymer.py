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
	"CH" : 1.09,
	"CO" : 1.41,
	"OC" : 1.41,
	"OCC" : 109.5,
	"COC" : 109.5,
	"CCO" : 109.5,
	"CCC" : 112.7
}

from utils import rad2deg,deg2rad,getTransformationMat
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

		# geometry and topology
		self.polyGRO_ = None
		self.polyITP_ = None

	def setNrepeats(self, Nrepeats):
		self.Nrepeats_ = Nrepeats

	def readPDB(self, PDBfile, preprocess=True):
		# read pdb file
		self.polyGRO_ = mda.Universe(PDBfile)

		if preprocess:
			# center the polymer
			self.center2graft()

			# align the minimum Rg axis to the x-axis (defaults)
			self.align(self.getOrient(), np.array([1, 0, 0]))

	def readGRO(self, GROfile):
		# read from files
		self.polyGRO_ = mda.Universe(GROfile)

		# center the polymer
		self.center2graft()
		# self.polyGRO_.atoms.write(f"test_centered.pdb")

		# align the minimum Rg axis to the x-axis (defaults)
		self.align(self.getOrient(), np.array([1, 0, 0]))
		# self.polyGRO_.atoms.write(f"test_aligned.pdb")

	def readITP(self, ITPfile):
		# read from files
		self.polyITP_ = mda.Universe(ITPfile, topology_format="ITP")

	def readDATA(self, DATAfile, atom_style="id resid type x y z"):
		# read lammps data file
		self.polyITP_ = mda.Universe(DATAfile, atom_style=atom_style)
		self.polyGRO_ = mda.Universe(DATAfile, atom_style=atom_style)

		# center the polymer
		self.center2graft()

		# align the minimum Rg axis to the x-axis (defaults)
		if self.polyname_ != "CG":
			self.align(self.getOrient(), np.array([1, 0, 0]))
		else:
			self.align(np.array([0, 0, 1]), np.array([1, 0, 0]))

		# set charge
		if "charge" not in atom_style.split():
			self.setCharge(0.0)

	def setCharge(self, charge=0.0):
		"""set the atom mass of the crystal"""
		atom_charges = [charge for _ in range(self.polyITP_.atoms.n_atoms)]
		self.polyGRO_.add_TopologyAttr('charge', atom_charges)
		self.polyITP_.add_TopologyAttr('charge', atom_charges)

	def getOrient(self):

		""" get the principal axis of a polymer chain for grafting use; 
		when polymer is aligned in a direction, 
		then the component in that direction is the smallest. 
		then the ref direction is the third principal axis, 
		the arrow is point from starting atom to the com"""
		
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
		# ref_vec = ref_vec/np.linalg.norm(ref_vec)
		return ref_vec

	def updatePos(self, pos):
		assert pos.shape[0] == self.polyGRO_.atoms.n_atoms, f"The given position should have the same length as the NO. of atoms!"
		assert pos.shape[1] == 3, f"The input position array should be three dimensional (x/y/z)! Exiting..."
		self.polyGRO_.load_new(pos, order='fac')

	def center2graft(self):
		"""
		move the grafting point (first atom) of a polymer to the origin 
		"""

		pos = self.polyGRO_.atoms.positions
		self.polyGRO_.load_new(pos - self.polyGRO_.atoms.positions[0], order='fac')

	def align(self, principal_axes, ref_direction):
		"""align one principal axis of a polymer to a given direction 
		https://docs.mdanalysis.org/1.0.0/documentation_pages/core/topologyattrs.html?highlight=principal%20axes#MDAnalysis.core.topologyattrs.Masses.principal_axes
		"""

		# # 0 is the longest axis
		# assert principal_axes_idx in [0,1,2], f"Only 0/1/2 of the axis is allowed for align the polymer!"
		# self.polyGRO_.atoms.align_principal_axis(principal_axes_idx, ref_direction)

		# use self-coded way
		RotMat,_ = getTransformationMat(principal_axes, ref_direction)
		updated_post = np.transpose(np.matmul(RotMat, np.transpose(np.array(self.polyGRO_.atoms.positions))))
		self.polyGRO_.load_new(updated_post, order='fac')

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

						# Hydrogen-s
						pos.append([xc,yc,zc+pars_dict['CH']])
						pos.append([xc,yc,zc-pars_dict['CH']])
						atomnames.append('H')
						atomnames.append('H')

						# second carbon
						xc = x_inc+dx
						yc = -dy
						zc = 0.0
						pos.append([xc,yc,zc])
						atomnames.append('C')

						# Hydrogen-s
						pos.append([xc,yc,zc+pars_dict['CH']])
						pos.append([xc,yc,zc-pars_dict['CH']])
						atomnames.append('H')
						atomnames.append('H')

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

						# Hydrogen-s
						pos.append([xc,yc,zc+pars_dict['CH']])
						pos.append([xc,yc,zc-pars_dict['CH']])
						atomnames.append('H')
						atomnames.append('H')

						# second carbon
						xc = x_inc+dx
						yc = 0.0
						zc = 0.0
						pos.append([xc,yc,zc])
						atomnames.append('C')

						# Hydrogen-s
						pos.append([xc,yc,zc+pars_dict['CH']])
						pos.append([xc,yc,zc-pars_dict['CH']])
						atomnames.append('H')
						atomnames.append('H')

						# oxygen
						xc = x_inc+2*dx
						yc = -dy
						zc = 0.0
						pos.append([xc,yc,zc])
						atomnames.append('O')

				# add last H with end O
				pos.append([xc,yc,zc+pars_dict['CH']])
				atomnames.append('H')

			else:
				print(f"Cyclic PEO is not available! Exiting...")
				sys.exit(0)

		elif self.polyname_ == "PVA":

			if topology == 'linear':
				# add head group
				xc = 0.0
				yc = 0.0
				zc = 0.0
				pos.append([xc, yc, zc])
				atomnames.append('C')

				# Hydrogen-s
				pos.append([xc,yc,zc+pars_dict['CH']])
				pos.append([xc,yc,zc-pars_dict['CH']])
				pos.append([xc-pars_dict['CH'],yc,zc])
				atomnames.append('H')
				atomnames.append('H')
				atomnames.append('H')

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

					# Hydrogen-s
					pos.append([xc,yc,zc+pars_dict['CH']])
					pos.append([xc,yc,zc-pars_dict['CH']])
					atomnames.append('H')
					atomnames.append('H')

					# second carbon
					xc = x_inc+2*dx
					yc = 0.0
					zc = 0.0
					pos.append([xc,yc,zc])
					atomnames.append('C')

					# Hydrogen-s
					pos.append([xc,yc,zc-pars_dict['CH']])
					atomnames.append('H')

					# oxygen
					pos.append([xc,yc+pars_dict["OC"],zc])
					atomnames.append('O')

					# Hydrogen-s
					pos.append([xc,yc+pars_dict["OC"],zc-pars_dict['CH']])
					atomnames.append('H')

				# add tail
				pos.append([xc+dx, -dy, 0.0])
				atomnames.append('C')

				# Hydrogen-s
				pos.append([xc+dx, -dy,pars_dict['CH']])
				pos.append([xc+dx, -dy,-pars_dict['CH']])
				pos.append([xc+dx+pars_dict['CH'],-dy,0.0])
				atomnames.append('H')
				atomnames.append('H')
				atomnames.append('H')

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

					# Hydrogen-s
					pos.append([xc,yc,zc+pars_dict['CH']])
					pos.append([xc,yc,zc-pars_dict['CH']])
					atomnames.append('H')
					atomnames.append('H')

					# second carbon
					xc = rad*math.cos((2*iatom+1)*theta)
					yc = rad*math.sin((2*iatom+1)*theta)
					zc = 0.0
					pos.append([xc,yc,zc])
					atomnames.append('C')

					# Hydrogen-s
					pos.append([xc,yc,zc-pars_dict['CH']])
					atomnames.append('H')

					# oxygen
					pos.append([xc,yc,zc+pars_dict["OC"]])
					atomnames.append('O')

					# Hydrogen-s
					pos.append([xc,yc,zc+pars_dict["OC"]+pars_dict["CH"]])
					atomnames.append('H')

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


class cgPolymer():

	def __init__(self):
		self.Natoms_ = 0
		self.atoms_ = None
		self.bonds_ = None
		self.angles_ = None
		self.dihedrals_ = None
		self.impropers_ = None

		self.n_atypes_ = 0
		self.n_btypes_ = 0
		self.positions_ = None

	def setChain(self, chain_length, block_frac, type_idx_shift=0):
		# block_frac is a list of fraction of each block
		assert np.sum(np.array(block_frac)) == 1, f"The total fraction should be 1!"

		self.Natoms_ = chain_length
		self.block_frac_ = block_frac
		self.n_atypes_ = len(block_frac)

		self.getTopo(type_idx_shift)
		self.getPos()

	def getAtomType(self):
		# get the type from the atom index
		block_idx = [0]
		for i in range(len(self.block_frac_)-1):
			start = sum(self.block_frac_[:i+1])*self.Natoms_
			block_idx.append(int(start))

		# add ends
		block_idx.append(self.Natoms_)

		# the out
		type_idx_all = [0 for _ in range(self.Natoms_)]
		atype = 0
		for iblock in range(len(block_idx)-1):
			atype += 1
			start = block_idx[iblock]
			end = block_idx[iblock+1]
			type_idx_all[start: end] = [atype for i in range(start, end)]

		return type_idx_all

	def bondidx2type(self, ibond):
		# get the bond type from the atom type of composing atoms

		a1 = ibond[0]
		a2 = ibond[1]

		# total atom types
		n_atypes = len(self.block_frac_)
		bond_dict = {}
		n_btypes = 0
		for atype_i in range(n_atypes):
			bond_dict[f"({atype_i+1}, {atype_i+1})"] = atype_i+1
			n_btypes += 1

		# cross term
		for atype_i in range(n_atypes-1):
			n_btypes += 1
			bond_dict[f"({atype_i+1}, {atype_i+2})"] = n_btypes

		self.n_btypes_ = len(bond_dict.keys())

		return bond_dict[f"({a1}, {a2})"]

	def getTopo(self, type_idx_shift=0):
		atomtype = [i+1 for i in range(len(self.block_frac_))]
		self.atoms_ = []
		self.bonds_ = []

		# atoms
		atype_idx_all = self.getAtomType()
		for iatom in range(self.Natoms_):
			atom_i = [iatom+1, 1, atype_idx_all[iatom]+type_idx_shift, 0, 0, iatom*1]
			self.atoms_.append(atom_i)

		# bonds
		for ibond in range(self.Natoms_-1):
			bond_i = [ibond+1, self.bondidx2type([atype_idx_all[ibond], atype_idx_all[ibond+1]]), ibond+1, ibond+2]
			self.bonds_.append(bond_i)

	def getPos(self):
		xyz = []

		for iatom in self.atoms_:
			xyz.append([iatom[3],iatom[4],iatom[5]])

		self.positions_ = np.array(xyz)

	def toDATA(self, fname):
		with open(fname, 'w') as FO:
			FO.write(f"lammps data file written by polyGraft\n")
			FO.write(f"\n")			
			FO.write(f"{self.Natoms_} atoms\n")
			FO.write(f"{self.Natoms_-1} bonds\n")
			FO.write(f"0 angles\n")			
			FO.write(f"0 dihedrals\n")			
			FO.write(f"0 impropers\n")			
			FO.write(f"\n")			
			FO.write(f"{self.n_atypes_} atom types\n")			
			FO.write(f"{self.n_btypes_} bond types\n")			
			FO.write(f"0 angle types\n")			
			FO.write(f"0 dihedral types\n")			
			FO.write(f"0 improper types\n")			

			pos = self.positions_
			box_max,box_min = np.amax(pos, axis=0), np.amin(pos, axis=0)
			FO.write(f"\n")
			FO.write(f"{box_min[0]-3:.3f} {box_max[0]+3:.3f} xlo xhi\n")
			FO.write(f"{box_min[1]-3:.3f} {box_max[1]+3:.3f} ylo yhi\n")
			FO.write(f"{box_min[2]-3:.3f} {box_max[2]+3:.3f} zlo zhi\n")

			# Atoms
			FO.write(f"\n")
			FO.write(f"Atoms\n")
			FO.write(f"\n")
			for iatom in self.atoms_:
				FO.write(f"{iatom[0]} {iatom[1]} {iatom[2]} {iatom[3]} {iatom[4]} {iatom[5]}\n")

			# Bonds
			FO.write(f"\n")
			FO.write(f"Bonds\n")
			FO.write(f"\n")
			for ibond in self.bonds_:
				FO.write(f"{ibond[0]} {ibond[1]} {ibond[2]} {ibond[3]}\n")

# if __name__ == '__main__':

	# # linear chains
	# linear = cgPolymer()
	# linear.setChain(12, [1.0], type_idx_shift=1)
	# linear.toDATA(f"linear_N12_type2.data")	
	
	# # block co-polymer
	# AB_bcp = cgPolymer()
	# AB_bcp.setChain(20, [0.5, 0.5])
	# AB_bcp.toDATA("bcp20_AB.data")

	# ABC_bcp = cgPolymer()
	# ABC_bcp.setChain(20, [0.3, 0.4, 0.3])
	# ABC_bcp.toDATA("bcp20_ABC.data")		