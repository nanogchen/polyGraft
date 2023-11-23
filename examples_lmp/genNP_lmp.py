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
	# lattice = Atomsk(lattice_type='fcc', lattice_const=4.08, element='Au')
	# lattice.gen_particle(radius, outFile="AuNP-R20.data")
	radius = 20.0
	nanoparticle = Crystal("nanoparticle", 'Au', radius)
	nanoparticle.readDATA("AuNP-R20.data",atom_style="id type x y z", guessing_bond=True, lattice_const=4.08)

	# graft
	peo_g_np = polyGraft(nanoparticle, peo)

	# set grafting density unit in A^-2
	gft = 0.0030 
	peo_g_np.setGraftingDensity(gft)

	# generate the grafted structure
	peo_g_np.setGftAtoms('Au')
	peo_g_np.genGraftStruct()

	# save gro and itp
	peo_g_np.toDATA("peo_g_np_gft"+"-R-"+str(radius)+"-sigma-"+str(gft)+".data")
