#!/bin/env python

import sys
sys.path.insert(0,"../polyGraft/")
from polyGraft import polyGraft
from polymer import Polymer
from cgCrystal import Crystal
from cgAtomsk import Atomsk

if __name__ == '__main__':

	# import linear
	linear = Polymer(poly_name="CG")

	# read from file
	chainLen = 12
	linear.readDATA(f"linear_N{chainLen}.data", atom_style="id resid type x y z")

	# read Au nanopore or by generation (the following two lines)
	lattice = Atomsk(lattice_type='fcc', nearest_neighbor=1, element='Au')
	radius = 5 # unit in lj (sigma)
	depth = 8
	lattice.gen_pore(radius, depth, outFile=f"Aupore.data")
	nanopore = Crystal("nanopore", 'Au', radius, depth)
	nanopore.readDATA(f"Aupore.data",atom_style="id type x y z", guessing_bond=True, nearest_neighbor=1)

	# graft
	poly_g_pore = polyGraft(nanopore, linear)

	# set grafting density unit in \sigma^-2
	gft_density = 0.1
	poly_g_pore.setGraftingDensity(gft_density)

	# generate the grafted structure
	poly_g_pore.setGftAtoms('Au')
	poly_g_pore.genGraftStruct()

	# save data
	poly_g_pore.toDATA(f"poly_g_pore_N{chainLen}-sigma-{str(gft_density)}-R{radius}.data", with_charges=False)

	
