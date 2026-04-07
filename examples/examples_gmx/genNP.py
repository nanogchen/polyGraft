#!/bin/env python

import sys
sys.path.insert(0,"../polyGraft/")
from polyGraft import polyGraft
from polymer import Polymer
from crystal import Crystal

if __name__ == '__main__':	

	# import peo
	peo = Polymer(poly_name="PEO")

	# read from file
	peo.readGRO("PEO12_line.gro")
	peo.readITP("PEO12.itp")

	# import Au nanoparticle
	radius = 20.0
	nanoparticle = Crystal("nanoparticle", 'Au', radius)
	nanoparticle.readPDB("AuNP-R20.pdb", guessing_bond=True, lattice_const=4.08)

	# graft
	peo_g_np = polyGraft(nanoparticle, peo)

	# set grafting density unit in A^-2
	gft = 0.0250 
	peo_g_np.setGraftingDensity(gft)

	# generate the grafted structure
	peo_g_np.setGftAtoms('Au')
	peo_g_np.genGraftStruct()

	# save gro and itp
	peo_g_np.toGRO("peo_g_np_gft"+"-R-"+str(radius)+"-sigma-"+str(gft)+".gro")
	peo_g_np.toITP("peo_g_np_gft"+"-R-"+str(radius)+"-sigma-"+str(gft)+".itp")
