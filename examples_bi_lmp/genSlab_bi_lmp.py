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
	peo.readDATA("../examples_lmp/PEO20.data", atom_style="id resid type charge x y z")

	# second graft
	peo6 = Polymer("PEO")
	peo6.readDATA("PEO6.data", atom_style="id resid type charge x y z")

	# read Au slab or by generation
	lattice = Atomsk(lattice_type='fcc', lattice_const=4.08, element='Au')
	length = 50.0
	width = 50.0
	depth = 10.0
	lattice.gen_slab(length,width,depth,outFile="Auslab.data")
	nanoslab = Crystal("nanoslab", 'Au', length, width, depth)
	nanoslab.readDATA("Auslab.data",atom_style="id type x y z", guessing_bond=True, lattice_const=4.08)

	# graft
	peo_g_slab = polyGraft(nanoslab, peo)
	peo_g_slab.setBinaryGraft(peo6)
	graft_type = 'homo-bigraft' 
	# graft_type = 'random-bigraft' 
	# graft_type = 'janus-bigraft' 
	peo_g_slab.setBinaryGraftStyle(graft_type)

	# set grafting density unit in A^-2
	gft = 0.0046 
	peo_g_slab.setGraftingDensity(gft)

	# generate the grafted structure
	peo_g_slab.setGftAtoms('Au')
	peo_g_slab.genGraftStruct()

	# save gro and itp
	peo_g_slab.toDATA(f"peo_g_slab_gft-sigma-{str(gft)}-{graft_type}_bi.data")
