#!/bin/env python

import sys
sys.path.insert(0, "../polyGraft/")
from polyGraft import polyGraft
from polymer import Polymer
from crystal import Crystal
from atomsk import Atomsk

if __name__ == '__main__':	

	# import peo
	peo = Polymer("PEO")

	# read
	peo.readGRO("PEO12_curve.gro")
	peo.readITP("PEO12.itp")

	# define lattice
	lattice = Atomsk(lattice_type='fcc', lattice_const=4.08, element='Au')
	
	# generate a pore
	radius = 30
	depth = 80
	lattice.gen_pore(radius, depth, outFile="Aupore.pdb")
	nanopore = Crystal("nanopore", 'Au', radius, depth)
	nanopore.readPDB("Aupore.pdb", guessing_bond=True, lattice_const=4.08)

	# graft
	peo_g_np = polyGraft(nanopore, peo)

	# set grafting density unit in A^-2
	gft = 0.0067
	peo_g_np.setGraftingDensity(gft)

	# generate the grafted structure
	peo_g_np.genGraftStruct()

	# save gro and itp
	peo_g_np.toGRO("peo_g_pore_gft"+"-R-"+str(radius)+"-sigma-"+str(gft)+".gro")
	peo_g_np.toITP("peo_g_pore_gft"+"-R-"+str(radius)+"-sigma-"+str(gft)+".itp")
