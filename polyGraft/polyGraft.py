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

import math
import sys
import numpy as np
import MDAnalysis as mda
from polymer import Polymer,pars_dict
from crystal import Crystal
from atomsk import Atomsk
from utils import getNN_two, getTransformationMat

class polyGraft():

	def __init__(self, center, graft):
		
		# member: center for grafts to graft
		self.center_ = center
		self.graft_ = graft
		self.graftingDensity_ = 0.0
		self.spacing_ = 0.0
		self.Ngrafts_ = 0
		self.Natoms_ = 0
		self.graftStruct_ = None
		self.graftAtoms_ = None # substrate atoms to be grafted with grafts
		self.graftPts_ = None
		self.centerIdx_ = None

	def setGraftingDensity(self, graftingDensity):
		if isinstance(self.center_, Crystal):
			assert graftingDensity < 0.06, f"Grafting density cannot be too large (smaller than 0.06 A^-2)!"
			self.graftingDensity_ = graftingDensity		

			# spacing 
			self.spacing_ = math.sqrt(1/graftingDensity)

		elif isinstance(self.center_, Polymer):
			assert self.graftingDensity_<=1.0, f"Grafting density ({self.graftingDensity_}) should not be larger than 1.0 grafts/monomer for bottlebrush polymers!"
			self.graftingDensity_ = graftingDensity

			# spacing 
			self.spacing_ = int(1.0/graftingDensity)

	def setGftAtoms(self, atomname):
		# set atoms of the substrate to be grafted
		if isinstance(self.center_, Crystal):
			self.graftAtoms_ = self.center_.crystal_.select_atoms(f'name {atomname}')

		elif isinstance(self.center_, Polymer):
			self.graftAtoms_ = self.center_.polyGRO_.select_atoms(f'name {atomname}')

	def genGraftStruct(self):
		# assemble two components together
		
		if isinstance(self.center_, Crystal):
			self.graftSoft2Hard()

		elif isinstance(self.center_, Polymer):
			self.graftSoft2Soft()

		else:
			print(f"Unexpected object type of the center object ({self.center_}): must be one of Polymer or Crystal.")
			sys.exit(0)

		# # merge grafts and sub
		# self.graftStruct_ = self.graft_.polyITP_

		# for igraft in range(1, self.Ngrafts_):
		# 	self.graftStruct_ = mda.Merge(self.graftStruct_.atoms, self.graft_.polyITP_.atoms)

		# self.graftStruct_ = mda.Merge(self.graftStruct_.atoms, self.center_.crystal_.atoms)

		# # set resid
		# resid_lst = [1 for _ in range(self.graftStruct_.atoms.n_atoms)]
		# self.graftStruct_.add_TopologyAttr('resid', resid_lst)

	def graftSoft2Soft(self):
		# for bottlebrush: generate the positions

		self.centerIdx_ = []
		self.graftStruct_ = self.center_.polyGRO_.copy()
		pos_all = []
		atomnames_all = []

		# loop over the graft points
		for i in range(self.graftAtoms_.n_atoms):

			# add side chain if satisfy the criteria
			if i % self.spacing_ == 0:
				self.centerIdx_.append(i+1)
				gft_pt = self.graftAtoms_[i].position
				i_gft = self.graft_.polyGRO_.copy()

				# two branches with 45 deg
				if i % 2 == 0:
					norm_vector = np.array([0, math.sqrt(2)/2.0, math.sqrt(2)/2.0])
				else:
					norm_vector = np.array([0, math.sqrt(2)/2.0, -math.sqrt(2)/2.0])

				ref_vect = self.graft_.getEndVect()
				RotMat,_ = getTransformationMat(ref_vect, norm_vector)
				TransMat = norm_vector*pars_dict['CC']
				pos_new = np.transpose(np.matmul(RotMat, np.transpose(np.array(i_gft.atoms.positions)))) + np.tile(gft_pt, ((len(i_gft.atoms.positions),1))) + np.tile(TransMat, ((len(i_gft.atoms.positions),1)))
				i_gft.load_new(pos_new, order='fac')

				self.graftStruct_ = mda.Merge(self.graftStruct_.atoms, i_gft.atoms)

	def graftSoft2Hard(self):
		
		if self.center_.crystal_type_ == "nanoparticle":
			self.particleGraft()			

		elif self.center_.crystal_type_ == "nanopore":
			self.poreGraft()

		elif self.center_.crystal_type_ == "nanorod":
			self.rodGraft()

		elif self.center_.crystal_type_ == "nanoslab":
			self.slabGraft()

		else:
			print(f"Unexpected crystal type of the center object ({self.center_.crystal_type_}): must be one of nanoparticle/nanopore/nanorod/nanoslab.")
			sys.exit(0)
	
	def slabGraft(self):
		# graft on the top surface of the slab

		length = self.center_.pars_[0]
		width = self.center_.pars_[1]
		depth = self.center_.pars_[2]

		# distribute chains in an unwrap rectangle
		Nchainsx = int(length/self.spacing_)
		Nchainsy = int(width/self.spacing_)
		self.Ngrafts_ = Nchainsx*Nchainsy

		# initialize the pos
		self.graftPts_ = np.zeros((self.Ngrafts_*self.graft_.polyGRO_.atoms.n_atoms,3))

		# self.Natoms
		self.Natoms_ = self.graft_.polyGRO_.atoms.n_atoms*self.Ngrafts_ + self.center_.crystal_.atoms.n_atoms

		# graft on the uppper side of the slab
		z = max(self.center_.crystal_.atoms.positions[:,2])
		x = [0.5*self.spacing_ + i*self.spacing_ for i in range(Nchainsx)]
		y = [0.5*self.spacing_ + i*self.spacing_ for i in range(Nchainsy)]
	
		graft_pts = []
		for i in x:
			for j in y:
				graft_pts.append([i,j,z])
		graft_pts = np.array(graft_pts)

		# find grafted AU pos/index
		gld_gft_pts, idx = getNN_two(graft_pts, self.center_.crystal_.atoms.positions)
		self.centerIdx_ = idx

		# determine normal and set
		poly_pos = self.graft_.polyGRO_.atoms.positions
		# ref_vect = self.graft_.getOrient()
		ref_vect = self.graft_.getEndVect()

		# get the pos of all chains
		for ichain in range(self.Ngrafts_):
			i_gft_pt =  gld_gft_pts[ichain,:]

			# get the normal vector
			norm_vector = np.array([0.0,0.0,1.0])

			RotMat,TransMat = getTransformationMat(ref_vect, norm_vector)
			RotMat = np.transpose(RotMat)

			ipos = np.matmul(poly_pos, RotMat) + np.tile(TransMat, (self.graft_.polyGRO_.atoms.n_atoms,1)) + np.tile(gld_gft_pts[ichain,:], ((self.graft_.polyGRO_.atoms.n_atoms,1)))
			self.graftPts_[ichain*self.graft_.polyGRO_.atoms.n_atoms:(1+ichain)*self.graft_.polyGRO_.atoms.n_atoms, :] = ipos

	def particleGraft(self):
		# graft on the outer surface of the particle

		radius = self.center_.pars_[0]

		# surface area
		surfArea = 4*math.pi*(radius**2)

		# calculate Nchains to graft
		self.Ngrafts_ = round(self.graftingDensity_*surfArea)

		# initialize the pos
		self.graftPts_ = np.zeros((self.Ngrafts_*self.graft_.polyGRO_.atoms.n_atoms,3))

		# self.Natoms
		self.Natoms_ = self.graft_.polyGRO_.atoms.n_atoms*self.Ngrafts_ + self.center_.crystal_.atoms.n_atoms

		# get graft pts
		center_pt = self.center_.crystal_.atoms.center_of_geometry()
		
		indices = np.arange(0, self.Ngrafts_, dtype=float) + 0.5
		phi = np.arccos(1 - 2*indices/self.Ngrafts_)
		theta = 2*np.pi * indices / (0.5*(1 + 5**0.5))
		x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
		
		# grafting pts: center at origin and move to the current
		graft_pts = np.array([x,y,z])*radius
		graft_pts = np.transpose(graft_pts) + np.expand_dims(np.array(center_pt), axis=0)

		# find gfted AU pos/index
		gld_gft_pts, idx = getNN_two(graft_pts, self.center_.crystal_.atoms.positions)
		self.centerIdx_ = idx

		# determine normal and set
		poly_pos = self.graft_.polyGRO_.atoms.positions
		# ref_vect = self.graft_.getOrient()
		ref_vect = self.graft_.getEndVect()

		# get the pos of all chains
		for ichain in range(self.Ngrafts_):
			i_gft_pt =  gld_gft_pts[ichain,:]

			# get the normal vector
			norm_vector = i_gft_pt - center_pt
			norm_vector = norm_vector/np.linalg.norm(norm_vector)

			RotMat,TransMat = getTransformationMat(ref_vect, norm_vector)
			RotMat = np.transpose(RotMat)

			ipos = np.matmul(poly_pos, RotMat) + np.tile(TransMat, (self.graft_.polyGRO_.atoms.n_atoms,1)) + np.tile(gld_gft_pts[ichain,:], ((self.graft_.polyGRO_.atoms.n_atoms,1)))
			self.graftPts_[ichain*self.graft_.polyGRO_.atoms.n_atoms:(1+ichain)*self.graft_.polyGRO_.atoms.n_atoms, :] = ipos

	def poreGraft(self):
		# graft on the inner surface of the pore

		radius = self.center_.pars_[0]
		depth = self.center_.pars_[1]

		# distribute chains in an unwrap rectangle
		length = 2*math.pi*radius
		Nchainsxy = int(length/self.spacing_)
		Nchainsz = int(depth/self.spacing_)
		self.Ngrafts_ = Nchainsxy*Nchainsz

		# initialize the pos
		self.graftPts_ = np.zeros((self.Ngrafts_*self.graft_.polyGRO_.atoms.n_atoms,3))

		# self.Natoms
		self.Natoms_ = self.graft_.polyGRO_.atoms.n_atoms*self.Ngrafts_ + self.center_.crystal_.atoms.n_atoms

		# graft on the inner surface of the pore
		center_pt = self.center_.crystal_.atoms.center_of_geometry()
		x = np.array([radius*math.cos(2*math.pi/Nchainsxy*i) for i in range(Nchainsxy)])
		y = np.array([radius*math.sin(2*math.pi/Nchainsxy*i) for i in range(Nchainsxy)])
		z = [0.5*self.spacing_ + i*self.spacing_ for i in range(Nchainsz)]

		graft_pts = []
		for i,j in zip(x,y):
			for k in z:
				graft_pts.append([i,j,k])

		# from local to global
		graft_pts = np.array(graft_pts)	+ center_pt - np.array([0,0,0.5*depth])

		# find grafted AU pos/index
		gld_gft_pts, idx = getNN_two(graft_pts, self.center_.crystal_.atoms.positions)
		self.centerIdx_ = idx

		# determine normal and set
		poly_pos = self.graft_.polyGRO_.atoms.positions
		# ref_vect = self.graft_.getOrient()
		# ref_vect = self.graft_.getEndVect()
		ref_vect = np.array([1.0,0.0,0.0])

		# get the pos of all chains
		for ichain in range(self.Ngrafts_):
			i_gft_pt =  gld_gft_pts[ichain,:]

			# get the normal vector
			norm_vector = center_pt - i_gft_pt
			norm_vector[2] = 0.0 # zero vertical componet
			norm_vector = norm_vector/np.linalg.norm(norm_vector)

			RotMat,TransMat = getTransformationMat(ref_vect, norm_vector)
			RotMat = np.transpose(RotMat)

			ipos = np.matmul(poly_pos, RotMat) + np.tile(TransMat, (self.graft_.polyGRO_.atoms.n_atoms,1)) + np.tile(gld_gft_pts[ichain,:], ((self.graft_.polyGRO_.atoms.n_atoms,1)))
			self.graftPts_[ichain*self.graft_.polyGRO_.atoms.n_atoms:(1+ichain)*self.graft_.polyGRO_.atoms.n_atoms, :] = ipos

	def rodGraft(self):
		# graft on the outer surface of the rod

		radius = self.center_.pars_[0]
		depth = self.center_.pars_[1]

		# distribute chains in an unwrap rectangle
		length = 2*math.pi*radius
		Nchainsxy = int(length/self.spacing_)
		Nchainsz = int(depth/self.spacing_)
		self.Ngrafts_ = Nchainsxy*Nchainsz

		# initialize the pos
		self.graftPts_ = np.zeros((self.Ngrafts_*self.graft_.polyGRO_.atoms.n_atoms,3))

		# self.Natoms
		self.Natoms_ = self.graft_.polyGRO_.atoms.n_atoms*self.Ngrafts_ + self.center_.crystal_.atoms.n_atoms

		# graft on the outer surface of the rod
		center_pt = self.center_.crystal_.atoms.center_of_geometry()
		x = np.array([radius*math.cos(2*math.pi/Nchainsxy*i) for i in range(Nchainsxy)])
		y = np.array([radius*math.sin(2*math.pi/Nchainsxy*i) for i in range(Nchainsxy)])
		z = [0.5*self.spacing_ + i*self.spacing_ for i in range(Nchainsz)]

		graft_pts = []
		for i,j in zip(x,y):
			for k in z:
				graft_pts.append([i,j,k])

		# from local to global
		graft_pts = np.array(graft_pts)	+ center_pt - np.array([0,0,0.5*depth])

		# find grafted AU pos/index
		gld_gft_pts, idx = getNN_two(graft_pts, self.center_.crystal_.atoms.positions)
		self.centerIdx_ = idx

		# determine normal and set
		poly_pos = self.graft_.polyGRO_.atoms.positions
		# ref_vect = self.graft_.getOrient()
		ref_vect = self.graft_.getEndVect()

		# get the pos of all chains
		for ichain in range(self.Ngrafts_):
			i_gft_pt =  gld_gft_pts[ichain,:]

			# get the normal vector
			norm_vector = i_gft_pt - center_pt
			norm_vector[2] = 0.0 # zero vertical componet
			norm_vector = norm_vector/np.linalg.norm(norm_vector)

			RotMat,TransMat = getTransformationMat(ref_vect, norm_vector)
			RotMat = np.transpose(RotMat)

			ipos = np.matmul(poly_pos, RotMat) + np.tile(TransMat, (self.graft_.polyGRO_.atoms.n_atoms,1)) + np.tile(gld_gft_pts[ichain,:], ((self.graft_.polyGRO_.atoms.n_atoms,1)))
			self.graftPts_[ichain*self.graft_.polyGRO_.atoms.n_atoms:(1+ichain)*self.graft_.polyGRO_.atoms.n_atoms, :] = ipos

	def toGRO(self, fname):

		""" save the assembled system in a gro file"""

		# write out gro file for VMD visualization
		with open(fname, 'w') as FO:

			# total atoms
			total_atoms = self.Natoms_

			# add header
			FO.write("Write by polyGraft \n")
			FO.write("%5d\n" % (total_atoms))

			# for iatom in self.graftStruct_.atoms:
			# 	FO.write(f"{iatom.resid:5d}{iatom.resname:5s}{iatom.name:5s}{iatom.id:5d}{iatom.position[0]/10:8.3f}{iatom.position[1]/10:8.3f}{iatom.position[2]/10:8.3f}\n")

			# write polymer first
			for ichain in range(self.Ngrafts_):
				for idx in range(self.graft_.polyGRO_.atoms.n_atoms):
					g_idx = ichain*self.graft_.polyGRO_.atoms.n_atoms + idx

					resid = self.graft_.polyGRO_.atoms.resids[idx]
					resname = self.graft_.polyGRO_.atoms.resnames[idx]
					atomname = self.graft_.polyGRO_.atoms.names[idx]
					atomid = g_idx + 1

					FO.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" % (resid, resname, atomname, atomid, \
								self.graftPts_[g_idx,0]/10, self.graftPts_[g_idx,1]/10, self.graftPts_[g_idx,2]/10))

			# write lattice 
			for idx in range(self.center_.crystal_.atoms.n_atoms):
				g_idx = self.Ngrafts_*self.graft_.polyGRO_.atoms.n_atoms + idx

				resid = idx + 1
				resname = self.center_.crystal_.atoms.resnames[idx]
				atomname = self.center_.crystal_.atoms.names[idx]
				atomid = g_idx + 1

				FO.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" % (resid, resname, atomname, atomid, \
							self.center_.crystal_.atoms.positions[idx][0]/10, \
							self.center_.crystal_.atoms.positions[idx][1]/10, \
							self.center_.crystal_.atoms.positions[idx][2]/10))

			# add ends
			FO.write("%10.5f%10.5f%10.5f\n" % (10.0,10.0,10.0))

	def toITP(self, fname):

		with open(fname, 'w') as FO:

			# header
			FO.write(f"[ moleculetype ]\n")
			FO.write(f"polyGraft            3\n")			
			FO.write("\n")

			# polymer grafts first, substrate second
			# atoms info
			FO.write("[ atoms ]\n")
			poly_n_atoms = self.graft_.polyITP_.atoms.n_atoms
			for i in range(self.Ngrafts_):
				for j in range(poly_n_atoms):
					atomidx = i*poly_n_atoms+j+1
					iatom = self.graft_.polyITP_.atoms[j]

					# atom_idx atom_type resid resname atom_name resid charge mass
					FO.write(f"{atomidx} {iatom.type} {iatom.resid} {iatom.resname} {iatom.name} {iatom.resid} {iatom.charge:.3f} {iatom.mass}\n")
			
			FO.write(";sub\n")
			for katom in self.center_.crystal_.atoms:
				atomidx+=1
				FO.write(f"{atomidx} {katom.type} {katom.resid} {katom.resname} {katom.name} {katom.resid} {katom.charge:.3f} {katom.mass}\n")
	
			# bonds info
			FO.write("\n")
			FO.write("[ bonds ]\n")
			for i in range(self.Ngrafts_):
				bond_shift = i*self.graft_.polyITP_.atoms.n_atoms

				for jbond in self.graft_.polyITP_.bonds:
					FO.write(f"{jbond._ix[0]+1+bond_shift} {jbond._ix[1]+1+bond_shift} {jbond.type}\n")
				
			FO.write(";sub-sub\n")
			bond_shift += self.graft_.polyITP_.atoms.n_atoms
			for kbond in self.center_.crystal_.bonds:
				# FO.write(f"{kbond._ix[0]+1+bond_shift} {kbond._ix[1]+1+bond_shift} {kbond.type}\n")
				FO.write(f"{kbond._ix[0]+1+bond_shift} {kbond._ix[1]+1+bond_shift} {1}\n")

			FO.write(";sub-graft\n")
			for ipoly in range(self.Ngrafts_):
				poly_idx = 1+ipoly*self.graft_.polyITP_.atoms.n_atoms
				FO.write(f"{poly_idx} {self.centerIdx_[ipoly]} {1}\n")

			# angles info
			FO.write("\n")
			FO.write("[ angles ]\n")
			for i in range(self.Ngrafts_):
				angle_shift = i*self.graft_.polyITP_.atoms.n_atoms
				
				for jangle in self.graft_.polyITP_.angles:					
					FO.write(f"{jangle._ix[0]+1+angle_shift} {jangle._ix[1]+1+angle_shift} {jangle._ix[2]+1+angle_shift} {jangle.type}\n")	

			# dihedrals info
			FO.write("\n")
			FO.write("[ dihedrals ]\n")
			for i in range(self.Ngrafts_):
				dihedral_shift = i*self.graft_.polyITP_.atoms.n_atoms

				for jdihedral in self.graft_.polyITP_.dihedrals:
					FO.write(f"{jdihedral._ix[0]+1+dihedral_shift} {jdihedral._ix[1]+1+dihedral_shift} {jdihedral._ix[2]+1+dihedral_shift} {jdihedral._ix[3]+1+dihedral_shift} {jdihedral.type}\n")
			
			# pairs info
			FO.write("\n")
			FO.write("[ pairs ]\n")
			FO.write(";from 1-4 pair\n")
			for i in range(self.Ngrafts_):
				dihedral_shift = i*self.graft_.polyITP_.atoms.n_atoms

				for jdihedral in self.graft_.polyITP_.dihedrals:
					FO.write(f"{jdihedral._ix[0]+1+dihedral_shift} {jdihedral._ix[3]+1+dihedral_shift} {1}\n")

	def toPDB(self, fname):

		if self.graftStruct_ is not None:
			self.graftStruct_.atoms.write(fname)

		else:
			print(f"The graft structure is None! Cannot write pdb file. Exiting")
			sys.exit(0)
