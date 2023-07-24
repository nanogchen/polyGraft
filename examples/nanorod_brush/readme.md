The generation example of nanorod brush. MUST use under /path/to/polyGraft/examples directory!

How to use:
```python
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
	peo.readGRO("PEO12_line.gro")
	peo.readITP("PEO12.itp")

	# define lattice
	lattice = Atomsk(lattice_type='fcc', lattice_const=4.08, element='Au')
	
	# generate a rod
	radius = 10.0
	depth = 50
	lattice.gen_rod(radius, depth, outFile="Aurod.pdb")
	nanorod = Crystal("nanorod", 'Au', radius, depth)
	nanorod.readPDB("Aurod.pdb", guessing_bond=True, lattice_const=4.08)

	# graft
	peo_g_np = polyGraft(nanorod, peo)

	# set grafting density unit in A^-2
	gft = 0.0250
	peo_g_np.setGraftingDensity(gft)

	# generate the grafted structure
	peo_g_np.setGftAtoms('Au')
	peo_g_np.genGraftStruct()

	# save gro and itp
	peo_g_np.toGRO("peo_g_rod_gft"+"-R-"+str(radius)+"-sigma-"+str(gft)+".gro")
	peo_g_np.toITP("peo_g_rod_gft"+"-R-"+str(radius)+"-sigma-"+str(gft)+".itp")
```