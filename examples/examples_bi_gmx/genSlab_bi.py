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
	peo.readGRO("../examples/PEO12_line.gro")
	peo.readITP("../examples/PEO12.itp")

	# second graft
	peo6 = Polymer("PEO")
	peo6.readGRO("PEO6.gro")
	peo6.readITP("PEO6.itp")

	# define lattice
	lattice = Atomsk(lattice_type='fcc', lattice_const=4.08, element='Au')

	# generate a slab
	length = 50.0
	width = 50.0
	depth = 10.0
	lattice.gen_slab(length,width,depth,outFile="Auslab.pdb")
	nanoslab = Crystal("nanoslab", 'Au', length, width, depth)
	nanoslab.readPDB("Auslab.pdb", guessing_bond=True, lattice_const=4.08)	

	# graft
	peo_g_np = polyGraft(nanoslab, peo)
	peo_g_np.setBinaryGraft(peo6)
	# graft_type = 'homo-bigraft' 
	graft_type = 'random-bigraft' 
	# graft_type = 'janus-bigraft' 
	peo_g_np.setBinaryGraftStyle(graft_type)

	# set grafting density unit in A^-2
	gft = 0.0120 
	peo_g_np.setGraftingDensity(gft)

	# generate the grafted structure
	peo_g_np.setGftAtoms('Au')
	peo_g_np.genGraftStruct()

	# save gro and itp
	peo_g_np.toGRO(f"peo_g_slab_gft-sigma-{str(gft)}-{graft_type}_bi.gro")
	peo_g_np.toITP(f"peo_g_slab_gft-sigma-{str(gft)}-{graft_type}_bi.itp")
