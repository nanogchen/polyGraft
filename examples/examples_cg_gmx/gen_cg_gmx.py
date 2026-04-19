#!/bin/env python

import sys,os
sys.path.insert(0, "../../src/")
from polyGraft import polyGraft
from polymer import Polymer
from crystal import Crystal
from atomsk import Atomsk

def gen(gro_file, itp_file, geometry, geom_params, lattice_type, lattice_constant, grafting_density, output_dir):
	
	# Debugging check:
	if not os.path.exists(gro_file):
		raise FileNotFoundError(f"Whoops! Expected to find {gro_file} but it isn't there.")
		
	if not os.path.exists(itp_file):
		raise FileNotFoundError(f"Whoops! Expected to find {itp_file} but it isn't there.")
		
	# import peo
	peo = Polymer("PEO")

	# read  
	peo.readGRO(gro_file)
	peo.readITP(itp_file)

	# define lattice
	lattice = Atomsk(lattice_type=lattice_type, lattice_const=lattice_constant, element='Au')

	# generate substrate
	tempfile = os.path.join(output_dir, "Au.pdb")	
	if geometry == "Slab":
		length = geom_params['Lx']
		width = geom_params['Ly']
		depth = geom_params['Lz']
		lattice.gen_slab(length,width,depth,outFile=tempfile)
		substrate = Crystal("nanoslab", 'Au', length, width, depth)
		substrate.readPDB(tempfile, guessing_bond=True, lattice_const=lattice_constant)

	elif geometry == "Rod":
		radius = geom_params['R']
		depth = geom_params['L']
		lattice.gen_rod(radius, depth, outFile=tempfile)
		substrate = Crystal("nanorod", 'Au', radius, depth)
		substrate.readPDB(tempfile, guessing_bond=True, lattice_const=lattice_constant)

	elif geometry == "Pore":
		radius = geom_params['R_in']
		depth = geom_params['L']
		lattice.gen_pore(radius, depth, outFile=tempfile)
		substrate = Crystal("nanopore", 'Au', radius, depth)
		substrate.readPDB(tempfile, guessing_bond=True, lattice_const=lattice_constant)

	elif geometry == "Sphere":
		radius = geom_params['Radius']
		lattice.gen_particle(radius, outFile=tempfile)
		substrate = Crystal("nanoparticle", 'Au', radius)
		substrate.readPDB(tempfile, guessing_bond=True, lattice_const=lattice_constant)	

	# graft
	peo_g_subs = polyGraft(substrate, peo)

	# set grafting density unit in A^-2
	peo_g_subs.setGraftingDensity(grafting_density)

	# generate the grafted structure
	peo_g_subs.setGftAtoms('Au')
	peo_g_subs.genGraftStruct()

	# Save the output to the provided temp output_dir
	final_gro_path = os.path.join(output_dir, "polygraft.gro")
	final_itp_path = os.path.join(output_dir, "polygraft.itp")
	peo_g_subs.toGRO(final_gro_path)
	peo_g_subs.toITP(final_itp_path)	
	
	return final_gro_path,final_itp_path 

if __name__ == '__main__':
	if len(sys.argv) != 9:
		print(f"Num of parameters is not correct!")
		sys.exit(0)

	gro_file = sys.argv[1]
	itp_file = sys.argv[2]
	geometry = sys.argv[3]
	geom_params = sys.argv[4]
	lattice_type = sys.argv[5]
	lattice_constant = float(sys.argv[6])
	grafting_density = float(sys.argv[7])
	output_dir = sys.argv[8]

	gen(gro_file, itp_file, geometry, geom_params, lattice_type, lattice_constant, grafting_density, output_dir)