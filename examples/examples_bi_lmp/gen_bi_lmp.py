#!/bin/env python

import sys,os
sys.path.insert(0, "../../src/")
from polyGraft import polyGraft
from polymer import Polymer
from crystal import Crystal
from atomsk import Atomsk

def gen(data1_file, data2_file, atom_style, geometry, geom_params, lattice_type, lattice_constant, bigraft_pattern, grafting_density, output_dir):
	
	# Debugging check:
	for datafile in [data1_file, data2_file]:
		if not os.path.exists(datafile):
			raise FileNotFoundError(f"Whoops! Expected to find {datafile} but it isn't there.")
		
	# import peo
	peo = Polymer("PEO")

	# read  
	peo.readDATA(data1_file, atom_style=atom_style)

	# second graft
	peo6 = Polymer("PEO")
	peo6.readDATA(data2_file, atom_style=atom_style)

	# define lattice
	lattice = Atomsk(lattice_type=lattice_type, lattice_const=lattice_constant, element='Au')

	# generate substrate
	tempfile = os.path.join(output_dir, "Au.data")	
	if geometry == "Slab":
		length = geom_params['Lx']
		width = geom_params['Ly']
		depth = geom_params['Lz']
		lattice.gen_slab(length,width,depth,outFile=tempfile)
		substrate = Crystal("nanoslab", 'Au', length, width, depth)
		substrate.readDATA(tempfile, atom_style="id type x y z", guessing_bond=True, lattice_const=lattice_constant)

	elif geometry == "Rod":
		radius = geom_params['R']
		depth = geom_params['L']
		lattice.gen_rod(radius, depth, outFile=tempfile)
		substrate = Crystal("nanorod", 'Au', radius, depth)
		substrate.readDATA(tempfile, atom_style="id type x y z", guessing_bond=True, lattice_const=lattice_constant)

	elif geometry == "Pore":
		radius = geom_params['R_in']
		depth = geom_params['L']
		lattice.gen_pore(radius, depth, outFile=tempfile)
		substrate = Crystal("nanopore", 'Au', radius, depth)
		substrate.readDATA(tempfile, atom_style="id type x y z", guessing_bond=True, lattice_const=lattice_constant)

	elif geometry == "Sphere":
		radius = geom_params['Radius']
		lattice.gen_particle(radius, outFile=tempfile)
		substrate = Crystal("nanoparticle", 'Au', radius)
		substrate.readDATA(tempfile, atom_style="id type x y z", guessing_bond=True, lattice_const=lattice_constant)	

	# graft
	peo_g_subs = polyGraft(substrate, peo)
	peo_g_subs.setBinaryGraft(peo6)
	peo_g_subs.setBinaryGraftStyle(bigraft_pattern)

	# set grafting density unit in A^-2
	peo_g_subs.setGraftingDensity(grafting_density)

	# generate the grafted structure
	peo_g_subs.setGftAtoms('Au')
	peo_g_subs.genGraftStruct()

	# Save the output to the provided temp output_dir
	final_data_path = os.path.join(output_dir, "polygraft.data")
	peo_g_subs.toDATA(final_data_path)
	
	return final_data_path,os.path.join(output_dir, "polygraft.prm")
