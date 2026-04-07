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
	peo.readDATA("PEO20.data", atom_style="id resid type charge x y z")

	# define lattice
	lattice = Atomsk(lattice_type='fcc', lattice_const=4.08, element='Au')
	
	# generate a rod
	radius = 10.0
	depth = 50
	lattice.gen_rod(radius, depth, outFile="Aurod.data")
	nanorod = Crystal("nanorod", 'Au', radius, depth)
	nanorod.readDATA("Aurod.data", guessing_bond=True, lattice_const=4.08)

	# graft
	peo_g_np = polyGraft(nanorod, peo)

	# set grafting density unit in A^-2
	gft = 0.0046
	peo_g_np.setGraftingDensity(gft)

	# generate the grafted structure
	peo_g_np.setGftAtoms('Au')
	peo_g_np.genGraftStruct()

	# save gro and itp
	peo_g_np.toDATA(f"peo_g_rod_gft-R-{radius}-sigma-{str(gft)}.data")
