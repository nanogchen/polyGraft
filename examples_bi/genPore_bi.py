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
	peo.readGRO("../examples/PEO12_curve.gro")
	peo.readITP("../examples/PEO12.itp")

	# second graft
	peo6 = Polymer("PEO")
	peo6.readGRO("PEO6.gro")
	peo6.readITP("PEO6.itp")

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
	peo_g_np.setBinaryGraft(peo6)
	# graft_type = 'homo-bigraft' 
	# graft_type = 'random-bigraft' 
	graft_type = 'janus-bigraft' 
	peo_g_np.setBinaryGraftStyle(graft_type)

	# set grafting density unit in A^-2
	gft = 0.0067 
	peo_g_np.setGraftingDensity(gft)

	# generate the grafted structure
	peo_g_np.setGftAtoms('Au')
	peo_g_np.genGraftStruct()

	# save gro and itp
	peo_g_np.toGRO(f"peo_g_pore_gft-R-{str(radius)}-sigma-{str(gft)}-{graft_type}_bi.gro")
	peo_g_np.toITP(f"peo_g_pore_gft-R-{str(radius)}-sigma-{str(gft)}-{graft_type}_bi.itp")
