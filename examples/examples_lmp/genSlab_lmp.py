#!/bin/env python

import sys
sys.path.insert(0,"../polyGraft/")
from polyGraft import polyGraft
from polymer import Polymer
from crystal import Crystal
from atomsk import Atomsk

if __name__ == '__main__':	

	# import peo
	peo = Polymer(poly_name="PEO")

	# read from file
	peo.readDATA("PEO20.data", atom_style="id resid type charge x y z")

	# read Au nanoparticle or by generation (the following two lines)
	lattice = Atomsk(lattice_type='fcc', lattice_const=4.08, element='Au')
	length = 50.0
	width = 50.0
	depth = 10.0
	lattice.gen_slab(length,width,depth,outFile="Auslab.data")
	nanoslab = Crystal("nanoslab", 'Au', length, width, depth)
	nanoslab.readDATA("Auslab.data",atom_style="id type x y z", guessing_bond=True, lattice_const=4.08)

	# graft
	peo_g_np = polyGraft(nanoslab, peo)

	# set grafting density unit in A^-2
	gft = 0.0046 
	peo_g_np.setGraftingDensity(gft)

	# generate the grafted structure
	peo_g_np.setGftAtoms('Au')
	peo_g_np.genGraftStruct()

	# save gro and itp
	peo_g_np.toDATA(f"peo_g_slab_gft-sigma-{str(gft)}.data")
