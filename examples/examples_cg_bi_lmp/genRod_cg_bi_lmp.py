#!/bin/env python

import sys
sys.path.insert(0,"../polyGraft/")
from polyGraft import polyGraft
from polymer import Polymer
from cgCrystal import Crystal
from cgAtomsk import Atomsk

if __name__ == '__main__':	

	# import first CG polymer
	poly1 = Polymer(poly_name="CG")
	poly1.readDATA("../examples_cg_lmp/linear_N12.data", atom_style="id resid type x y z")

	# second graft
	poly2 = Polymer("CG")
	poly2.readDATA("linear_N6_shifted.data", atom_style="id resid type x y z")

	# read Au nanorod or by generation (the following two lines)
	lattice = Atomsk(lattice_type='fcc', nearest_neighbor=1, element='Au')
	radius = 4 # unit in lj (sigma)
	depth = 15
	lattice.gen_rod(radius, depth, outFile=f"Aurod.data")
	nanorod = Crystal("nanorod", 'Au', radius, depth)
	nanorod.readDATA(f"Aurod.data",atom_style="id type x y z", guessing_bond=True, nearest_neighbor=1)

	# graft
	poly_g_rod = polyGraft(nanorod, poly1)
	poly_g_rod.setBinaryGraft(poly2)
	# graft_type = 'homo-bigraft' 
	# graft_type = 'random-bigraft' 
	graft_type = 'janus-bigraft' 
	poly_g_rod.setBinaryGraftStyle(graft_type)

	# set grafting density unit in A^-2
	gft = 0.20 
	poly_g_rod.setGraftingDensity(gft)

	# generate the grafted structure
	poly_g_rod.setGftAtoms('Au')
	poly_g_rod.genGraftStruct()

	# save gro and itp
	poly_g_rod.toDATA(f"poly_g_rod_gft-sigma-{str(gft)}-{graft_type}_bi.data", with_charges=False)
