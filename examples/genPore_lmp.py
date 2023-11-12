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
	
	# generate a pore
	radius = 50
	depth = 40
	lattice.gen_pore(radius, depth, outFile="Aupore.data")
	nanopore = Crystal("nanopore", 'Au', radius, depth)
	nanopore.readDATA("Aupore.data", guessing_bond=True, lattice_const=4.08)

	# graft
	peo_g_np = polyGraft(nanopore, peo)

	# set grafting density unit in A^-2
	gft = 0.0030
	peo_g_np.setGraftingDensity(gft)

	# generate the grafted structure
	peo_g_np.setGftAtoms('Au')
	peo_g_np.genGraftStruct()

	# save gro and itp3
	peo_g_np.toDATA(f"peo_g_pore_gft-R-{radius}-sigma-{str(gft)}.data")
