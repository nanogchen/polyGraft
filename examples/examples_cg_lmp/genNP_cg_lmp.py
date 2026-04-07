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

	# read Au nanoparticle or by generation (the following two lines)
	lattice = Atomsk(lattice_type='fcc', nearest_neighbor=1, element='Au')
	radius = 5 # unit in lj (sigma)
	lattice.gen_particle(radius,outFile=f"AuNP{radius}.data")
	nanoparticle = Crystal("nanoparticle", 'Au', radius)
	nanoparticle.readDATA(f"AuNP{radius}.data",atom_style="id type x y z", guessing_bond=True, nearest_neighbor=1)

	# graft
	poly_g_np = polyGraft(nanoparticle, linear)

	# set grafting density unit in \sigma^-2
	gft_density = 0.1
	poly_g_np.setGraftingDensity(gft_density)

	# generate the grafted structure
	poly_g_np.setGftAtoms('Au')
	poly_g_np.genGraftStruct()

	# save data
	poly_g_np.toDATA(f"poly_g_np_N12-sigma-{str(gft_density)}-R{radius}.data", with_charges=False)

	
