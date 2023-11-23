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

import random
import math
import sys
sys.path.insert(0, "../examples/")
from rtp_define import gen_BBP_rtp, res_rtp_dict
import numpy as np
import MDAnalysis as mda
from polymer import Polymer,pars_dict
from crystal import Crystal
from atomsk import Atomsk
from utils import getNN_two, getTransformationMat
import utils

class polyGraft():

	def __init__(self, center, graft):
		
		# member: center for grafts to graft
		self.center_ = center
		self.graft_ = graft
		self.graftingDensity_ = 0.0
		self.spacing_ = 0.0
		self.Ngrafts_ = 0
		self.graftStruct_ = None	# the final hybrid structure
		self.graftAtoms_ = None 	# substrate atoms optional for grafts

		self.all_grafts_lst_ = None # list of all grafts (atom group)
		self.graft_chain_idx_ = None# list of chain idx for each graft point
		self.centerGftedIdx_ = None # substrate atoms to be grafted with grafts

		# atom groups for poly-g-soft/BBPs
		self.centerAtmGrp_ = None
		self.graftAtmGrp_ = None

		# for bigrafts
		self.bigraft_mode_  = None

	def setBinaryGraft(self, graft_bi):
		# read in the second graft for binary grafts
		self.graft_bi_ = graft_bi
		
	def setBinaryGraftStyle(self, style):
		# set the grafting style for binary grafts
		self.bigraft_mode_ = style
		assert style in ['homo-bigraft','random-bigraft','janus-bigraft'], f"Currently only homo/random/janus type is allowed!"

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

			# check the grafting density
			assert self.center_.Nrepeats_ % self.spacing_ == 0, f"For bottlebrush generation, the backbone length ({self.center_.Nrepeats_}) should be N times of the spacing ({self.spacing_})!"

	def setGftAtoms(self, atomname):
		# set atoms of the substrate to be grafted
		if isinstance(self.center_, Crystal):
			self.graftAtoms_ = self.center_.crystal_.select_atoms(f'name {atomname}')

			# the final structure poly-g-soft: initialized with the backbones
			self.graftStruct_ = self.center_.crystal_.atoms
			self.centerAtmGrp_ = self.center_.crystal_.atoms

		elif isinstance(self.center_, Polymer):
			self.graftAtoms_ = self.center_.polyGRO_.select_atoms(f'name {atomname}')

			# find the index of backbone atoms with grafts
			self.centerGftedIdx_ = []
			for i in range(self.graftAtoms_.n_atoms):
				if i % self.spacing_ == 0:
					self.centerGftedIdx_.append(self.graftAtoms_[i].index)

			# the total number of grafts
			self.Ngrafts_ = len(self.centerGftedIdx_)

			# remove Hydrogen bonded to graft atom if there is a graft
			index_lst = list(set(self.center_.polyGRO_.atoms.indices)-set(np.array(self.centerGftedIdx_)+1))
			
			# the final structure poly-g-soft: initialized with the backbones
			self.graftStruct_ = self.center_.polyGRO_.atoms[index_lst]
			self.centerAtmGrp_ = self.center_.polyGRO_.atoms[index_lst]
			
	def genGraftStruct(self, random_grafting=False):
		# assemble two components together
		
		if isinstance(self.center_, Crystal):
			self.graftSoft2Hard(random_grafting)

		elif isinstance(self.center_, Polymer):
			self.graftSoft2Soft(random_grafting)

		else:
			print(f"Unexpected object type of the center object ({self.center_}): must be one of Polymer or Crystal.")
			sys.exit(0)

	def graftSoft2Soft(self, random_grafting=False):
		# for bottlebrush: generate the positions

		pos_all = []
		atomnames_all = []

		# loop over the graft points
		ngrafts = 0
		for i in range(self.graftAtoms_.n_atoms):

			# add side chain if satisfy the criteria
			if i % self.spacing_ == 0:
				ngrafts += 1
				gft_pt = self.graftAtoms_[i].position
				i_gft = self.graft_.polyGRO_.copy()

				# two branches with 45 deg
				if i % 2 == 0:
					if self.center_.topology_ == "linear":
						norm_vector = np.array([0, math.sqrt(2)/2.0, math.sqrt(2)/2.0])
					else:
						norm_vector = np.array([0, 0, 1])

				else:
					if self.center_.topology_ == "linear":
						norm_vector = np.array([0, math.sqrt(2)/2.0, -math.sqrt(2)/2.0])
					else:
						norm_vector = np.array([0, 0, 1])

				ref_vect = self.graft_.getEndVect()
				RotMat,_ = getTransformationMat(ref_vect, norm_vector)
				TransMat = norm_vector*pars_dict['CC']
				pos_new = np.transpose(np.matmul(RotMat, np.transpose(np.array(i_gft.atoms.positions)))) + np.tile(gft_pt, ((len(i_gft.atoms.positions),1))) + np.tile(TransMat, ((len(i_gft.atoms.positions),1)))
				i_gft.load_new(pos_new, order='fac')

				self.graftStruct_ = mda.Merge(self.graftStruct_.atoms, i_gft.atoms)

				if ngrafts == 1:
					self.graftAtmGrp_ = i_gft.copy()

				else:
					self.graftAtmGrp_ = mda.Merge(self.graftAtmGrp_.atoms, i_gft.atoms)

	def graftSoft2Hard(self, random_grafting):
		
		if self.center_.crystal_type_ == "nanoparticle":
			self.particleGraft()			

		elif self.center_.crystal_type_ == "nanopore":
			self.poreGraft()

		elif self.center_.crystal_type_ == "nanorod":
			self.rodGraft()

		elif self.center_.crystal_type_ == "nanoslab":
			self.slabGraft(random_grafting)

		else:
			print(f"Unexpected crystal type of the center object ({self.center_.crystal_type_}): must be one of nanoparticle/nanopore/nanorod/nanoslab.")
			sys.exit(0)
	
	def slabGraft(self, random_grafting=False):
		# graft on the top surface of the slab

		length = self.center_.pars_[0]
		width = self.center_.pars_[1]
		depth = self.center_.pars_[2]

		# distribute chains in an unwrap rectangle
		Nchainsx = int(length/self.spacing_)
		Nchainsy = int(width/self.spacing_)
		self.Ngrafts_ = Nchainsx*Nchainsy

		# initialize the pos
		self.all_grafts_lst_ = []

		# graft on the upper side of the slab
		z = max(self.center_.crystal_.atoms.positions[:,2])
		graft_pts = []

		if random_grafting:
			for _ in range(self.Ngrafts_):
				random_x = random.random()
				x = length*random_x
				random_y = random.random()
				y = width*random_y
				graft_pts.append([x,y,z])

		else:
			x = [0.5*self.spacing_ + i*self.spacing_ for i in range(Nchainsx)]
			y = [0.5*self.spacing_ + i*self.spacing_ for i in range(Nchainsy)]

			for i in x:
				for j in y:
					graft_pts.append([i,j,z])

		graft_pts = np.array(graft_pts)

		# find grafted AU pos/index
		gld_gft_pts, idx = getNN_two(graft_pts, self.graftAtoms_.positions)
		self.centerGftedIdx_ = idx

		# obtain idx [0 for 1st graft; 1 for 2nd graft]
		self.graft_chain_idx_ = [0 for _ in range(self.Ngrafts_)]
		if self.bigraft_mode_ == "homo-bigraft":
			for igft in range(self.Ngrafts_):
				if igft % 2 != 0:
					self.graft_chain_idx_[igft] = 1

		elif self.bigraft_mode_ == "random-bigraft":
			random_idx = random.choices([i for i in range(self.Ngrafts_)], k=int(self.Ngrafts_/2))
			for i in random_idx:
				self.graft_chain_idx_[i] = 1

		elif self.bigraft_mode_ == "janus-bigraft":
			for igft in range(self.Ngrafts_):
				if gld_gft_pts[igft, 0] > 0.5*length :
					self.graft_chain_idx_[igft] = 1

		# get the pos of all chains
		for ichain in range(self.Ngrafts_):
			i_gft_pt =  gld_gft_pts[ichain,:]

			# get the normal vector
			norm_vector = np.array([0.0,0.0,1.0])

			# select graft based on the selection (determine normal and set)
			if self.graft_chain_idx_[ichain] == 0:
				poly_pos = self.graft_.polyGRO_.atoms.positions
				ref_vect = self.graft_.getEndVect()
				i_graft = self.graft_.polyGRO_.copy()

			else:
				poly_pos = self.graft_bi_.polyGRO_.atoms.positions
				ref_vect = self.graft_bi_.getEndVect()
				i_graft = self.graft_bi_.polyGRO_.copy()

			RotMat,TransMat = getTransformationMat(ref_vect, norm_vector)
			RotMat = np.transpose(RotMat)

			ipos = np.matmul(poly_pos, RotMat) + np.tile(TransMat, (i_graft.atoms.n_atoms,1)) + np.tile(gld_gft_pts[ichain,:], ((i_graft.atoms.n_atoms,1)))
			i_graft.load_new(ipos, order='fac')
			self.all_grafts_lst_.append(i_graft)

			self.graftStruct_ = mda.Merge(self.graftStruct_.atoms, i_graft.atoms)

			if ichain == 0:
				self.graftAtmGrp_ = i_graft.copy()

			else:
				self.graftAtmGrp_ = mda.Merge(self.graftAtmGrp_.atoms, i_graft.atoms)

	def particleGraft(self):
		# graft on the outer surface of the particle

		radius = self.center_.pars_[0]

		# surface area
		surfArea = 4*math.pi*(radius**2)

		# calculate Nchains to graft
		self.Ngrafts_ = round(self.graftingDensity_*surfArea)

		# initialize the pos
		self.grafts_all_pos_ = np.zeros((self.Ngrafts_*self.graft_.polyGRO_.atoms.n_atoms,3))

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
		gld_gft_pts, idx = getNN_two(graft_pts, self.graftAtoms_.positions)
		self.centerGftedIdx_ = idx

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
			self.grafts_all_pos_[ichain*self.graft_.polyGRO_.atoms.n_atoms:(1+ichain)*self.graft_.polyGRO_.atoms.n_atoms, :] = ipos

			i_gft = self.graft_.polyGRO_.copy()
			i_gft.load_new(ipos, order='fac')

			self.graftStruct_ = mda.Merge(self.graftStruct_.atoms, i_gft.atoms)

			if ichain == 0:
				self.graftAtmGrp_ = i_gft.copy()

			else:
				self.graftAtmGrp_ = mda.Merge(self.graftAtmGrp_.atoms, i_gft.atoms)

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
		self.grafts_all_pos_ = np.zeros((self.Ngrafts_*self.graft_.polyGRO_.atoms.n_atoms,3))

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
		gld_gft_pts, idx = getNN_two(graft_pts, self.graftAtoms_.positions)
		self.centerGftedIdx_ = idx

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
			self.grafts_all_pos_[ichain*self.graft_.polyGRO_.atoms.n_atoms:(1+ichain)*self.graft_.polyGRO_.atoms.n_atoms, :] = ipos

			i_gft = self.graft_.polyGRO_.copy()
			i_gft.load_new(ipos, order='fac')

			self.graftStruct_ = mda.Merge(self.graftStruct_.atoms, i_gft.atoms)

			if ichain == 0:
				self.graftAtmGrp_ = i_gft.copy()

			else:
				self.graftAtmGrp_ = mda.Merge(self.graftAtmGrp_.atoms, i_gft.atoms)

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
		self.grafts_all_pos_ = np.zeros((self.Ngrafts_*self.graft_.polyGRO_.atoms.n_atoms,3))

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
		gld_gft_pts, idx = getNN_two(graft_pts, self.graftAtoms_.positions)
		self.centerGftedIdx_ = idx

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
			self.grafts_all_pos_[ichain*self.graft_.polyGRO_.atoms.n_atoms:(1+ichain)*self.graft_.polyGRO_.atoms.n_atoms, :] = ipos

			i_gft = self.graft_.polyGRO_.copy()
			i_gft.load_new(ipos, order='fac')

			self.graftStruct_ = mda.Merge(self.graftStruct_.atoms, i_gft.atoms)

			if ichain == 0:
				self.graftAtmGrp_ = i_gft.copy()

			else:
				self.graftAtmGrp_ = mda.Merge(self.graftAtmGrp_.atoms, i_gft.atoms)
				
	def toGRO(self, fname):

		""" save the assembled system in a gro file"""

		if isinstance(self.center_, Polymer):
			if self.graftStruct_ is not None:
				self.graftStruct_.atoms.write(fname)

			else:
				print(f"The graft structure is None! Cannot write gro file. Exiting")
				sys.exit(0)			

		else:			

			# write out gro file for VMD visualization
			with open(fname, 'w') as FO:

				# total atoms
				total_atoms = self.graftStruct_.atoms.n_atoms

				# add header
				FO.write("Write by polyGraft \n")
				FO.write("%5d\n" % (total_atoms))

				# write polymer first
				NgraftsAtoms = 0
				for igraft in self.all_grafts_lst_:
					ipos = igraft.atoms.positions
					for idx in range(igraft.atoms.n_atoms):

						resid = igraft.atoms.resids[idx]
						resname = igraft.atoms.resnames[idx]
						atomname = igraft.atoms.names[idx]
						atomid = NgraftsAtoms + idx + 1

						FO.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" % (resid, resname, atomname, atomid, \
									ipos[idx,0]/10, ipos[idx,1]/10, ipos[idx,2]/10))

					NgraftsAtoms += igraft.atoms.n_atoms

				# write lattice 
				for idx in range(self.center_.crystal_.atoms.n_atoms):
					g_idx = NgraftsAtoms + idx

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

		if isinstance(self.center_, Polymer):
			print(f"Generation of itp file for bottlebrush polymer is not available by this function! Generate pdb and rtp then use Gromacs to generate it!")
			sys.exit(0)			

		else:	

			with open(fname, 'w') as FO:

				# header
				FO.write(f"[ moleculetype ]\n")
				FO.write(f"polyGraft            3\n")			
				FO.write("\n")

				# polymer grafts first, substrate second
				# atoms info
				FO.write("[ atoms ]\n")
				NgraftsAtoms = 0
				for graft_chain in self.graft_chain_idx_:
					if graft_chain == 0:
						igraft = self.graft_
					else:
						igraft = self.graft_bi_

					poly_n_atoms = igraft.polyITP_.atoms.n_atoms
					for idx in range(poly_n_atoms):

						atomidx = NgraftsAtoms+idx+1
						iatom = igraft.polyITP_.atoms[idx]

						# atom_idx atom_type resid resname atom_name resid charge mass
						FO.write(f"{atomidx} {iatom.type} {iatom.resid} {iatom.resname} {iatom.name} {iatom.resid} {iatom.charge:.3f} {iatom.mass}\n")
					
					NgraftsAtoms += igraft.polyITP_.atoms.n_atoms

				FO.write(";sub\n")
				for katom in self.center_.crystal_.atoms:
					atomidx+=1
					FO.write(f"{atomidx} {katom.type} {katom.resid} {katom.resname} {katom.name} {katom.resid} {katom.charge:.3f} {katom.mass}\n")
		
				# bonds info
				FO.write("\n")
				FO.write("[ bonds ]\n")
				bond_shift = 0
				for graft_chain in self.graft_chain_idx_:
					if graft_chain == 0:
						igraft = self.graft_
					else:
						igraft = self.graft_bi_

					for jbond in igraft.polyITP_.bonds:
						FO.write(f"{jbond._ix[0]+1+bond_shift} {jbond._ix[1]+1+bond_shift} {jbond.type}\n")

					bond_shift += igraft.polyITP_.atoms.n_atoms
					
				FO.write(";sub-sub\n")
				for kbond in self.center_.crystal_.bonds:
					# FO.write(f"{kbond._ix[0]+1+bond_shift} {kbond._ix[1]+1+bond_shift} {kbond.type}\n")
					FO.write(f"{kbond._ix[0]+1+bond_shift} {kbond._ix[1]+1+bond_shift} {1}\n")

				FO.write(";sub-graft\n")
				atom_shift = 0
				for idx, graft_chain in enumerate(self.graft_chain_idx_):
					if graft_chain == 0:
						igraft = self.graft_
					else:
						igraft = self.graft_bi_

					poly_idx = 1+atom_shift
					FO.write(f"{poly_idx} {self.centerGftedIdx_[idx]+NgraftsAtoms} {1}\n")

					atom_shift += igraft.polyITP_.atoms.n_atoms

				# angles info: for now only consider polymer
				FO.write("\n")
				FO.write("[ angles ]\n")
				angle_shift = 0
				for graft_chain in self.graft_chain_idx_:
					if graft_chain == 0:
						igraft = self.graft_
					else:
						igraft = self.graft_bi_

					for jangle in igraft.polyITP_.angles:					
						FO.write(f"{jangle._ix[0]+1+angle_shift} {jangle._ix[1]+1+angle_shift} {jangle._ix[2]+1+angle_shift} {jangle.type}\n")	

					angle_shift += igraft.polyITP_.atoms.n_atoms

				# dihedrals info
				FO.write("\n")
				FO.write("[ dihedrals ]\n")
				dihedral_shift = 0
				for graft_chain in self.graft_chain_idx_:
					if graft_chain == 0:
						igraft = self.graft_
					else:
						igraft = self.graft_bi_	

					for jdihedral in igraft.polyITP_.dihedrals:
						FO.write(f"{jdihedral._ix[0]+1+dihedral_shift} {jdihedral._ix[1]+1+dihedral_shift} {jdihedral._ix[2]+1+dihedral_shift} {jdihedral._ix[3]+1+dihedral_shift} {jdihedral.type}\n")
					
					dihedral_shift += igraft.polyITP_.atoms.n_atoms

				# pairs info
				FO.write("\n")
				FO.write("[ pairs ]\n")
				FO.write(";from 1-4 pair\n")
				dihedral_shift = 0
				for graft_chain in self.graft_chain_idx_:
					if graft_chain == 0:
						igraft = self.graft_
					else:
						igraft = self.graft_bi_

					for jdihedral in self.graft_.polyITP_.dihedrals:
						FO.write(f"{jdihedral._ix[0]+1+dihedral_shift} {jdihedral._ix[3]+1+dihedral_shift} {1}\n")

					dihedral_shift += igraft.polyITP_.atoms.n_atoms

	def toDATA(self, fname):
		# write result to lammps data file

		# change the atom type in letter to integers
		atomtypes = self.graftStruct_.atoms.types
		unique_atypes = list(set(atomtypes))
		unique_atypes.sort()

		# consider graft-subs cross bond + 1
		unique_btypes = utils.getUniqueBondTypes(self.graftStruct_)
		unique_btypes.sort()
		# print(unique_btypes)

		# angles: no cross
		unique_angtypes = utils.getUniqueAngleTypes(self.graftStruct_)
		unique_angtypes.sort()

		# dihedrals: no cross
		unique_dtypes = utils.getUniqueDihedralTypes(self.graftStruct_)
		unique_dtypes.sort()

		# impropers: no cross
		unique_itypes = utils.getUniqueImproperTypes(self.graftStruct_)
		unique_itypes.sort()

		# define the atom types by idx
		atomtype2index = {k:v+1 for v,k in enumerate(unique_atypes)}

		# define the bond types by idx
		bondtype2index = {k:v+1 for v,k in enumerate(unique_btypes)}
		angletype2index = {k:v+1 for v,k in enumerate(unique_angtypes)}
		dihedraltype2index = {k:v+1 for v,k in enumerate(unique_dtypes)}
		impropertype2index = {k:v+1 for v,k in enumerate(unique_itypes)}

		# write the dictionary
		with open("polyGraft.prm", 'w') as fo:
			fo.write(f"# atoms\n")
			for key, value in atomtype2index.items():
				fo.write(f"{key}:{value}\n")

			fo.write(f"\n")
			fo.write(f"# bonds\n")
			for key, value in bondtype2index.items():
				fo.write(f"{key}:{value}\n")

			fo.write(f"\n")
			fo.write(f"# angles\n")
			for key, value in angletype2index.items():
				fo.write(f"{key}:{value}\n")

			fo.write(f"\n")
			fo.write(f"# dihedrals\n")
			for key, value in dihedraltype2index.items():
				fo.write(f"{key}:{value}\n")

			fo.write(f"\n")
			fo.write(f"# impropers\n")
			for key, value in impropertype2index.items():
				fo.write(f"{key}:{value}\n")

		# in-house write data out file
		with open(fname, 'w') as FO:

			# add header
			FO.write("lammps data file written by polyGraft \n")
			FO.write(f"\n")			
			FO.write(f"{self.graftStruct_.atoms.n_atoms} atoms\n")
			FO.write(f"{len(self.center_.crystal_.bonds)+self.Ngrafts_+self.Ngrafts_*len(self.graft_.polyGRO_.bonds)} bonds\n")
			FO.write(f"{len(self.graftStruct_.atoms.angles)} angles\n")
			FO.write(f"{len(self.graftStruct_.atoms.dihedrals)} dihedrals\n")
			FO.write(f"{len(self.graftStruct_.atoms.impropers)} impropers\n")
			FO.write(f"\n")
			FO.write(f"{len(unique_atypes)} atom types\n") 
			FO.write(f"{len(unique_btypes)+1} bond types\n")
			FO.write(f"{len(unique_angtypes)} angle types\n")
			FO.write(f"{len(unique_dtypes)} dihedral types\n")
			FO.write(f"{len(unique_itypes)} improper types\n")

			pos = self.graftStruct_.atoms.positions
			box_max,box_min = np.amax(pos, axis=0), np.amin(pos, axis=0)
			FO.write(f"\n")
			FO.write(f"{box_min[0]-3:.3f} {box_max[0]+3:.3f} xlo xhi\n")
			FO.write(f"{box_min[1]-3:.3f} {box_max[1]+3:.3f} ylo yhi\n")
			FO.write(f"{box_min[2]-3:.3f} {box_max[2]+3:.3f} zlo zhi\n")

			# Atoms
			FO.write(f"\n")
			FO.write(f"Atoms\n")
			FO.write(f"\n")

			# write polymer first
			atomtypes_vec = self.graft_.polyGRO_.atoms.types
			indexed_atomtypes_vec = [atomtype2index[i] for i in atomtypes_vec]
			NgraftsAtoms = 0
			for ichain, igraft in enumerate(self.all_grafts_lst_):
				ipos = igraft.atoms.positions
				for atom_idx in range(igraft.atoms.n_atoms):
					g_idx = NgraftsAtoms + atom_idx

					FO.write(f"{g_idx+1} {ichain+1} {indexed_atomtypes_vec[atom_idx]} {igraft.atoms.charges[atom_idx]:.3f} {ipos[atom_idx,0]:.3f} {ipos[atom_idx,1]:.3f} {ipos[atom_idx,2]:.3f}\n")
				
				NgraftsAtoms += igraft.atoms.n_atoms

			# write lattice 
			ichain+=1
			for idx in range(self.center_.crystal_.atoms.n_atoms):
				g_idx = NgraftsAtoms + idx

				FO.write(f"{g_idx+1} {ichain+1} {atomtype2index[self.center_.crystal_.atoms[0].type]} {self.center_.crystal_.atoms.charges[atom_idx]:.3f} {self.center_.crystal_.atoms.positions[idx,0]:.3f} {self.center_.crystal_.atoms.positions[idx,1]:.3f} {self.center_.crystal_.atoms.positions[idx,2]:.3f}\n")

			# Bonds
			FO.write(f"\n")
			FO.write(f"Bonds\n")
			FO.write(f"\n")
			# polymer first
			nbonds = 0
			bond_shift = 0
			for graft_chain in self.graft_chain_idx_:
				if graft_chain == 0:
					igraft = self.graft_
				else:
					igraft = self.graft_bi_

				for ibond in igraft.polyITP_.bonds:
					nbonds += 1
					FO.write(f"{nbonds} {ibond.type} {ibond._ix[0]+1+bond_shift} {ibond._ix[1]+1+bond_shift}\n")
				
				bond_shift += igraft.polyITP_.atoms.n_atoms

			# lattice
			for jbond in self.center_.crystal_.bonds:
				nbonds += 1
				FO.write(f"{nbonds} {bondtype2index[jbond.type]} {jbond._ix[0]+1+bond_shift} {jbond._ix[1]+1+bond_shift}\n")

			# polymer-lattice
			atom_shift = 0
			for idx, graft_chain in enumerate(self.graft_chain_idx_):
				if graft_chain == 0:
					igraft = self.graft_
				else:
					igraft = self.graft_bi_

				poly_idx = 1+atom_shift			
				nbonds += 1
				FO.write(f"{nbonds} {len(bondtype2index.keys())+1} {poly_idx} {self.centerGftedIdx_[idx]+NgraftsAtoms}\n")
				atom_shift += igraft.polyITP_.atoms.n_atoms

			# angles
			FO.write(f"\n")
			FO.write(f"Angles\n")
			FO.write(f"\n")
			nangles = 0
			angle_shift = 0
			for graft_chain in self.graft_chain_idx_:
				if graft_chain == 0:
					igraft = self.graft_
				else:
					igraft = self.graft_bi_
	
				for jangle in igraft.polyITP_.angles:
					nangles += 1
					FO.write(f"{nangles} {jangle.type} {jangle._ix[0]+1+angle_shift} {jangle._ix[1]+1+angle_shift} {jangle._ix[2]+1+angle_shift}\n")	

				angle_shift += igraft.polyITP_.atoms.n_atoms

			# dihedrals
			FO.write(f"\n")
			FO.write(f"Dihedrals\n")
			FO.write(f"\n")
			ndihedrals = 0
			dihedral_shift = 0
			for graft_chain in self.graft_chain_idx_:
				if graft_chain == 0:
					igraft = self.graft_
				else:
					igraft = self.graft_bi_	

				for jdihedral in igraft.polyITP_.dihedrals:
					ndihedrals += 1
					FO.write(f"{ndihedrals} {jdihedral.type} {jdihedral._ix[0]+1+dihedral_shift} {jdihedral._ix[1]+1+dihedral_shift} {jdihedral._ix[2]+1+dihedral_shift} {jdihedral._ix[3]+1+dihedral_shift}\n")
				
				dihedral_shift += igraft.polyITP_.atoms.n_atoms

			# impropers
			FO.write(f"\n")
			FO.write(f"Impropers\n")
			FO.write(f"\n")
			nimpropers = 0
			improper_shift = 0
			for graft_chain in self.graft_chain_idx_:
				if graft_chain == 0:
					igraft = self.graft_
				else:
					igraft = self.graft_bi_
	
				for jimproper in igraft.polyITP_.impropers:
					nimpropers += 1
					FO.write(f"{nimpropers} {jimproper.type} {jimproper._ix[0]+1+improper_shift} {jimproper._ix[1]+1+improper_shift} {jimproper._ix[2]+1+improper_shift} {jimproper._ix[3]+1+improper_shift}\n")

				improper_shift += igraft.polyITP_.atoms.n_atoms

	def toPDB(self, fname):

		if self.graftStruct_ is not None:

			# write out pdb file: wrap up atom/residue names
			self.PDBwrap(fname)

			# topology needed to add 
			self.checkTOP(self.center_.topology_)

		else:
			print(f"The graft structure is None! Cannot write pdb file. Exiting")
			sys.exit(0)

	def PDBwrap(self, fname):
	# wrap up pdb files of updated molecular file (after added Hydrogen in Avogadro)

		BTB_G_atoms, _ = gen_BBP_rtp(self.graft_.Nrepeats_, self.spacing_)

		# get res for 1/2/3
		BTB2_G = BTB_G_atoms.copy()
		if self.center_.topology_ == 'linear':
			BTB1_G = ['C31','HE1','HE2','HE3']
			BTB1_G.extend(BTB_G_atoms)
			BTB3_G = BTB_G_atoms.copy()
			BTB3_G.extend(['C31','HE1','HE2','HE3'])
		else:
			BTB1_G = BTB_G_atoms.copy()
			BTB3_G = BTB_G_atoms.copy()

		# read in the updated pdb file
		C_pos_bb = self.centerAtmGrp_.select_atoms('name C').positions.tolist()
		O_pos_bb = self.centerAtmGrp_.select_atoms('name O').positions.tolist()
		H_pos_bb = self.centerAtmGrp_.select_atoms('name H').positions.tolist()

		C_pos_sc = self.graftAtmGrp_.select_atoms('name C').positions.tolist()
		O_pos_sc = self.graftAtmGrp_.select_atoms('name O').positions.tolist()
		H_pos_sc = self.graftAtmGrp_.select_atoms('name H').positions.tolist()

		# bb or sc atoms
		bb_lst = ['C3', 'C2', 'CT','OH', 'HE', 'HT','HP','HO']
		sc_lst = ['CS', 'OP', 'OA', 'HA','H1','H2','H3','H4','H5','H6','H7','H8','H9'] 

		# use a dict to save the xyz
		all_pos_bb = {'O':O_pos_bb,\
				   'C':C_pos_bb,\
				   'H':H_pos_bb}

		all_pos_sc = {'O':O_pos_sc,\
				   'C':C_pos_sc,\
				   'H':H_pos_sc}

		# residues: BTB1 (1) + BTB2 (N-2) + BTB3 (1)	
		with open(fname, 'w') as FO:

			# header
			FO.write("TITLE Write by polyGraft-pdbwrap \n")

			# atomidx
			atomid = 0

			# out variables
			rectype = 'ATOM'

			# loop over residues
			for i_Gft in range(self.Ngrafts_):
				resid = i_Gft+1

				# BTB1
				if i_Gft == 0:

					for atom_idx in range(len(BTB1_G)):
						atomid += 1
						atomname =  BTB1_G[atom_idx]
						resname = 'BTB1'
						if atomname[:2] in bb_lst:
							xyz = all_pos_bb[BTB1_G[atom_idx][0]].pop(0)
						elif atomname[:2] in sc_lst:
							if atomname == 'OP1':
								xyz = all_pos_bb[BTB1_G[atom_idx][0]].pop(0)
							else:
								xyz = all_pos_sc[BTB1_G[atom_idx][0]].pop(0)

						occ = 	  1.00
						tempfac = 0.00
						element = BTB1_G[atom_idx][0]

						FO.write("{:<6s}{:>5d} {:>4s} {:>3s} {:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>12s}\n".format(\
								rectype, atomid, atomname, resname, resid, xyz[0], xyz[1], xyz[2],occ,tempfac,element))

				# BTB3
				elif i_Gft == self.Ngrafts_-1:

					for atom_idx in range(len(BTB3_G)):
						atomid += 1
						atomname =  BTB3_G[atom_idx]
						resname = 'BTB3'
						if atomname[:2] in bb_lst:
							xyz = all_pos_bb[BTB3_G[atom_idx][0]].pop(0)
						elif atomname[:2] in sc_lst:
							if atomname == 'OP1':
								xyz = all_pos_bb[BTB3_G[atom_idx][0]].pop(0)
							else:
								xyz = all_pos_sc[BTB3_G[atom_idx][0]].pop(0)

						occ = 	  1.00
						tempfac = 0.00
						element = BTB3_G[atom_idx][0]

						FO.write("{:<6s}{:>5d} {:>4s} {:>3s} {:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>12s}\n".format(\
								rectype, atomid, atomname, resname, resid, xyz[0], xyz[1], xyz[2],occ,tempfac,element))

				# BTB2
				else:

					for atom_idx in range(len(BTB2_G)):
						atomid += 1
						atomname =  BTB2_G[atom_idx]
						resname = 'BTB2'
						if atomname[:2] in bb_lst:
							xyz = all_pos_bb[BTB2_G[atom_idx][0]].pop(0)
						elif atomname[:2] in sc_lst:
							if atomname == 'OP1':
								xyz = all_pos_bb[BTB2_G[atom_idx][0]].pop(0)
							else:
								xyz = all_pos_sc[BTB2_G[atom_idx][0]].pop(0)

						occ = 	  1.00
						tempfac = 0.00
						element = BTB2_G[atom_idx][0]

						FO.write("{:<6s}{:>5d} {:>4s} {:>3s} {:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>12s}\n".format(\
								rectype, atomid, atomname, resname, resid, xyz[0], xyz[1], xyz[2],occ,tempfac,element))

			# end
			FO.write("END\n")

	def checkTOP(self, topology):
		# topology needed to add for cyclic bbps

		if topology == "cyclic":
			fname = "Nbb"+str(self.center_.Nrepeats_)+"_Nsc"+str(self.graft_.Nrepeats_)+"_sigma"+str(self.graftingDensity_)+"_to_add_cyclic.top"
			with open(fname, "w") as fo:

				# index of head and tail
				head = 1
				Natoms_per_bb_monomer = self.center_.polyGRO_.atoms.n_atoms/self.center_.Nrepeats_

				if self.centerGftedIdx_[-1] > (self.center_.Nrepeats_-1)*Natoms_per_bb_monomer:
					# last c grafted
					# backbone atoms - (Hs removed) + side chain atoms | atoms after the last monomer
					tail = Natoms_per_bb_monomer*(self.center_.Nrepeats_-1) - (self.Ngrafts_-1) + self.graft_.polyGRO_.atoms.n_atoms*(self.Ngrafts_-1) + 4 
				else:
					tail = Natoms_per_bb_monomer*(self.center_.Nrepeats_-1) - (self.Ngrafts_) + self.graft_.polyGRO_.atoms.n_atoms*self.Ngrafts_ + 4

				fo.write(f";bonds head-to-tail\n")
				fo.write(f"{head:>5d}{int(tail):>6d}{1:>6d}\n")
				fo.write(f"\n")
				fo.write(f";angles head-to-tail\n")
				fo.write(f"to add!\n")
				fo.write(f"\n")
				fo.write(f";dihedrals head-to-tail\n")
				fo.write(f"to add!\n")
				fo.write(f"\n")
				fo.write(f";pairs 1-4 of dihedrals\n")
				fo.write(f"to add!\n")

	def toRTP(self, fname):
		# generate the residue definition

		if isinstance(self.center_, Crystal):
			print(f"Generation of rtp file for poly-g-hard brush is not available by this function! Directly generate gro and itp file!")
			sys.exit(0)			

		else:
			BBP_G_atoms, BBP_G_bonds = gen_BBP_rtp(self.graft_.Nrepeats_, self.spacing_)

			with open(fname, 'w') as FO:
				# header
				FO.write("[ bondedtypes ]\n")
				FO.write("; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih\n")
				FO.write("     1       1          3          1	    1         3      1     0\n")
				FO.write("\n")

				for i in range(3):
					FO.write(f"[ BTB{i+1} ]\n")

					# ------------------------------ write atoms
					FO.write(" [ atoms ]\n")

					# BTB1
					if i==0:
						if self.center_.topology_ == "linear":
							FO.write(f"; CH3\n")
							FO.write(f"C3{1} {res_rtp_dict['C3']}\n")
							FO.write(f"HE{1} {res_rtp_dict['HE']}\n")
							FO.write(f"HE{2} {res_rtp_dict['HE']}\n")
							FO.write(f"HE{3} {res_rtp_dict['HE']}\n")

					# BTB2
					for iatom in BBP_G_atoms:
						try: # if not H
							FO.write(f"{iatom} {res_rtp_dict[iatom[:2]]}\n")
						except:
							FO.write(f"{iatom} {res_rtp_dict['H']}\n")

					# BTB3
					if i==2:
						if self.center_.topology_ == "linear":
							FO.write(f"; CH3\n")
							FO.write(f"C3{1} {res_rtp_dict['C3']}\n")
							FO.write(f"HE{1} {res_rtp_dict['HE']}\n")
							FO.write(f"HE{2} {res_rtp_dict['HE']}\n")
							FO.write(f"HE{3} {res_rtp_dict['HE']}\n")

					# ------------------------------ write bonds
					FO.write(" [ bonds ]\n")
					if i==0:
						if self.center_.topology_ == "linear":
							# add head
							FO.write(f"C3{1} {BBP_G_atoms[0]}\n")
							FO.write(f"C3{1} HE{1}\n")
							FO.write(f"C3{1} HE{2}\n")
							FO.write(f"C3{1} HE{3}\n")
							FO.write(f"{BBP_G_bonds[-1]}\n")

						else:
							# with next residue
							FO.write(f"{BBP_G_bonds[-1]}\n")

					if i==2:
						# with previous residue
						FO.write(f"{BBP_G_bonds[-2]}\n")

					for ibond in BBP_G_bonds[:-2]:
						FO.write(f"{ibond}\n")

					if i==1:
						FO.write(f"{BBP_G_bonds[-1]}\n")
						FO.write(f"{BBP_G_bonds[-2]}\n")

					if i==2:
						if self.center_.topology_ == "linear":
							# add head
							FO.write(f"C3{1} CT{self.spacing_}\n")
							FO.write(f"C3{1} HE{1}\n")
							FO.write(f"C3{1} HE{2}\n")
							FO.write(f"C3{1} HE{3}\n")				

					# add new line if the res is done
					FO.write("\n")
