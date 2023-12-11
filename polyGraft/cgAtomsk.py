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
import math
import os
import sys
sys.path.insert(0, '../bin/')

class Atomsk:

	def __init__(self, lattice_type, nearest_neighbor, element):
		assert lattice_type in ['fcc', 'bcc'], f"Currently only fcc/bcc lattice type is supported, but {lattice_type} was given!"

		self.lattice_type_ = lattice_type
		self.nearest_neighbor_ = nearest_neighbor		
		self.element_ = element

		if lattice_type == 'fcc':
			self.lattice_const_ = math.sqrt(2)
		elif lattice_type == 'bcc':
			self.lattice_const_ = 2/math.sqrt(3)

		# crystal type
		self.crystal_type_ = None

		# print information
		print(f"Please know that length unit in cgAtomsk is in 1-sigma of LJ unit!")
		assert os.path.exists("../bin/atomsk"), f"The atomsk program should be placed under /path/to/polyGraft/bin/ to generate atomic lattice crystals"

	def xsf2pdb(self, infname):
		# check if the file already exists
		pdbfile = infname.split(".")[0]+".pdb"
		self.clean_files([pdbfile])

		# convert xsf to pdb
		os.system(f"atomsk {infname} pdb >> gen.log")

		# clean
		self.clean_files(["gen.log"])

	def xsf2data(self, infname):
		# check if the file already exists
		pre = infname.split(".")[0]
		datafile = pre + ".data"
		self.clean_files([datafile])

		# convert xsf to pdb
		os.system(f"atomsk {infname} lammps >> gen.log")
		os.system(f"mv {pre}.lmp {pre}.data")

		# clean
		self.clean_files(["gen.log"])

	def pdb2gro(self, infname):
		# pdb cannot be directly read in MDAnalysis
		pre = infname.split(".")[0]
		outgro = pre+".gro"
		self.clean_files([outgro])

		# convert pdb to gro
		os.system(f"gmx editconf -f {pre}.pdb -o {outgro} > gmx.log")

		# clean
		self.clean_files([pre+".pdb", "gmx.log"])

	def clean_files(self, file_list):

		for ifile in file_list:
			if os.path.exists(ifile):
				os.remove(ifile)

	def create_unit_cell(self, fname):

		# check if the file already exists
		if os.path.exists(fname):
			os.remove(fname)

		# use atomsk
		os.system(f"atomsk --create {self.lattice_type_} {self.lattice_const_} {self.element_} {fname} ")
			
	def gen_slab(self, length=10.0, width=10.0, depth=4.0, outFile="slab.pdb"):
		# get the number of duplicates in each direction
		Nx = int(length/self.lattice_const_)
		Ny = int(width/self.lattice_const_)
		Nz = int(depth/self.lattice_const_)

		# check if the xsf file already exists
		xsffile = outFile.split(".")[0]+".xsf"
		self.clean_files([xsffile])

		# use atomsk to generate the lattice
		os.system(f"atomsk --create {self.lattice_type_} {self.lattice_const_} {self.element_} -duplicate {Nx} {Ny} {Nz} {xsffile} >> gen.log")

		# convert xsf to pdb
		outFileFormat = outFile.split(".")[-1]
		if outFileFormat == "pdb":
			self.xsf2pdb(xsffile)
		elif outFileFormat == "data":
			self.xsf2data(xsffile)
		else:
			print(f"Unknown out file format (outFileFormat) found! Can only be pdb or lammps-data...")
			sys.exit(0)

		# clean outFiles
		self.clean_files([xsffile, "gen.log"])

		# set type
		self.crystal_type_ = "slab"

	def gen_particle(self, radius, outFile='particle.data'):

		# get the number of duplicates
		Ndup = int(2*radius/self.lattice_const_)

		# check if the cubic xsf file already exists
		xsffile = "cubic.xsf"
		self.clean_files([xsffile])

		# generate a cubic lattice
		os.system(f"atomsk --create {self.lattice_type_} {self.lattice_const_} {self.element_} -duplicate {Ndup} {Ndup} {Ndup} {xsffile}  >> gen.log")

		# check if the cubic xsf file already exists
		outxsffile = outFile.split(".")[0] + ".xsf"
		self.clean_files([outxsffile])

		# cut the cubic to form particle
		os.system(f"atomsk {xsffile} -select out sphere 0.5*box 0.5*box 0.5*box {radius} -rmatom select {outxsffile}  >> gen.log")

		# convert xsf to pdb
		outFileFormat = outFile.split(".")[1]
		if outFileFormat == "pdb":
			self.xsf2pdb(outxsffile)
		elif outFileFormat == "data":
			self.xsf2data(outxsffile)
		else:
			print(f"Unknown out file format (outFileFormat) found! Can only be pdb or lammps-data...")
			sys.exit(0)

		# clean outFiles
		self.clean_files([xsffile, outxsffile, "gen.log"])

		# set type
		self.crystal_type_ = "particle"

	def gen_pore(self, radius, depth, outFile='pore.pdb'):

		# get the parameter for a slab first
		Nxy = int(2*radius/self.lattice_const_)+2
		Nz = int(depth/self.lattice_const_)

		# check if the cubic xsf file already exists
		xsffile = "cubic.xsf"
		self.clean_files([xsffile])

		# generate a cubic lattice
		os.system(f"atomsk --create {self.lattice_type_} {self.lattice_const_} {self.element_} -duplicate {Nxy} {Nxy} {Nz} {xsffile} >> gen.log")

		# check if the cubic xsf file already exists
		outxsffile = outFile.split(".")[0] + ".xsf"
		self.clean_files([outxsffile])

		# cut the cubic to form particle
		os.system(f"atomsk {xsffile} -select in cylinder Z 0.49*box 0.49*box {radius} -rmatom select {outxsffile} >> gen.log")

		# convert xsf to pdb
		outFileFormat = outFile.split(".")[1]
		if outFileFormat == "pdb":
			self.xsf2pdb(outxsffile)
		elif outFileFormat == "data":
			self.xsf2data(outxsffile)
		else:
			print(f"Unknown out file format (outFileFormat) found! Can only be pdb or lammps-data...")
			sys.exit(0)

		# clean outFiles
		self.clean_files([xsffile, outxsffile, "gen.log"])

		# set type
		self.crystal_type_ = "pore"

	def gen_rod(self, radius, depth, outFile='rod.pdb'):

		# get the parameter for a slab first
		Nxy = int(2*radius/self.lattice_const_)+3
		Nz = int(depth/self.lattice_const_)

		# check if the cubic xsf file already exists
		xsffile = "cubic.xsf"
		self.clean_files([xsffile])

		# generate a cubic lattice
		os.system(f"atomsk --create {self.lattice_type_} {self.lattice_const_} {self.element_} -duplicate {Nxy} {Nxy} {Nz} {xsffile} >> gen.log")

		# check if the cubic xsf file already exists
		outxsffile = outFile.split(".")[0] + ".xsf"
		self.clean_files([outxsffile])

		# cut the cubic to form particle
		os.system(f"atomsk {xsffile} -select out cylinder Z 0.5*box 0.5*box {radius} -rmatom select {outxsffile} >> gen.log")

		# convert xsf to pdb
		outFileFormat = outFile.split(".")[1]
		if outFileFormat == "pdb":
			self.xsf2pdb(outxsffile)
		elif outFileFormat == "data":
			self.xsf2data(outxsffile)
		else:
			print(f"Unknown out file format (outFileFormat) found! Can only be pdb or lammps-data...")
			sys.exit(0)

		# clean outFiles
		self.clean_files([xsffile, outxsffile, "gen.log"])

		# set type
		self.crystal_type_ = "rod"

