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
	linear.readDATA(f"linear_N12.data", atom_style="id resid type x y z")

	# read Au nanoslab or by generation (the following two lines)
	lattice = Atomsk(lattice_type='fcc', nearest_neighbor=1, element='Au')
	length = 10 # unit in lj (sigma)
	width = 10
	depth = 4
	lattice.gen_slab(length,width,depth,outFile="Auslab.data")
	nanoslab = Crystal("nanoslab", 'Au', length, width, depth)
	nanoslab.readDATA("Auslab.data",atom_style="id type x y z", guessing_bond=True, nearest_neighbor=1)

	# graft
	poly_g_slab = polyGraft(nanoslab, linear)

	# set grafting density unit in \sigma^-2
	gft_density = 0.1
	poly_g_slab.setGraftingDensity(gft_density)

	# generate the grafted structure
	poly_g_slab.setGftAtoms('Au')
	poly_g_slab.genGraftStruct()

	# save data
	poly_g_slab.toDATA(f"poly_g_slab_N12-sigma-{str(gft_density)}.data", with_charges=False)
	
