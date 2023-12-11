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
	linear.readDATA(f"linear_N6.data", atom_style="id resid type x y z")

	# read Au nanorod or by generation (the following two lines)
	lattice = Atomsk(lattice_type='fcc', nearest_neighbor=1, element='Au')
	radius = 4 # unit in lj (sigma)
	depth = 15
	lattice.gen_rod(radius, depth, outFile=f"Aurod.data")
	nanorod = Crystal("nanorod", 'Au', radius, depth)
	nanorod.readDATA(f"Aurod.data",atom_style="id type x y z", guessing_bond=True, nearest_neighbor=1)

	# graft
	poly_g_rod = polyGraft(nanorod, linear)

	# set grafting density unit in \sigma^-2
	gft_density = 0.1
	poly_g_rod.setGraftingDensity(gft_density)

	# generate the grafted structure
	poly_g_rod.setGftAtoms('Au')
	poly_g_rod.genGraftStruct()

	# save data
	poly_g_rod.toDATA(f"poly_g_rod_N6-sigma-{str(gft_density)}-R{radius}.data", with_charges=False)

	
