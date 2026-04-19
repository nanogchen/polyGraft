#!/bin/env python

import sys,os
sys.path.insert(0, "../../src/")
from polyGraft import polyGraft
from polymer import Polymer
from cgCrystal import Crystal
from cgAtomsk import Atomsk

def gen(data_file, atom_style, geometry, geom_params, lattice_type, lattice_constant, grafting_density, output_dir):
	
	# Debugging check:
	if not os.path.exists(data_file):
		raise FileNotFoundError(f"Whoops! Expected to find {gro_file} but it isn't there.")
		
	# import peo
	peo = Polymer("CG")

	# read  
	peo.readDATA(data_file, atom_style=atom_style)

	# define lattice
	lattice = Atomsk(lattice_type=lattice_type, nearest_neighbor=1, element='Au')

	# generate substrate
	tempfile = os.path.join(output_dir, "Au.data")	
	if geometry == "Slab":
		length = geom_params['Lx']
		width = geom_params['Ly']
		depth = geom_params['Lz']
		lattice.gen_slab(length,width,depth,outFile=tempfile)
		substrate = Crystal("nanoslab", 'Au', length, width, depth)
		substrate.readDATA(tempfile, atom_style="id type x y z", guessing_bond=True, nearest_neighbor=lattice_constant)

	elif geometry == "Rod":
		radius = geom_params['R']
		depth = geom_params['L']
		lattice.gen_rod(radius, depth, outFile=tempfile)
		substrate = Crystal("nanorod", 'Au', radius, depth)
		substrate.readDATA(tempfile, atom_style="id type x y z", guessing_bond=True, nearest_neighbor=lattice_constant)

	elif geometry == "Pore":
		radius = geom_params['R_in']
		depth = geom_params['L']
		lattice.gen_pore(radius, depth, outFile=tempfile)
		substrate = Crystal("nanopore", 'Au', radius, depth)
		substrate.readDATA(tempfile, atom_style="id type x y z", guessing_bond=True, nearest_neighbor=lattice_constant)

	elif geometry == "Sphere":
		radius = geom_params['Radius']
		lattice.gen_particle(radius, outFile=tempfile)
		substrate = Crystal("nanoparticle", 'Au', radius)
		substrate.readDATA(tempfile, atom_style="id type x y z", guessing_bond=True, nearest_neighbor=lattice_constant)	

	# graft
	peo_g_subs = polyGraft(substrate, peo)

	# set grafting density unit in A^-2
	peo_g_subs.setGraftingDensity(grafting_density)

	# generate the grafted structure
	peo_g_subs.setGftAtoms('Au')
	peo_g_subs.genGraftStruct()

	# Save the output to the provided temp output_dir
	final_data_path = os.path.join(output_dir, "polygraft.data")
	peo_g_subs.toDATA(final_data_path, with_charges=False)
	
	return final_data_path,os.path.join(output_dir, "polygraft.prm")
