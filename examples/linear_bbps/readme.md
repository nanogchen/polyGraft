The generation example of linear bottlebrush polymer. MUST use under /path/to/polyGraft/examples directory!

How to use:
```python
#!/bin/env python

import sys
sys.path.insert(0, "../polyGraft/")
from polymer import Polymer
from polyGraft import polyGraft

if __name__ == '__main__':	

	# generate PEO side chain
	peo=Polymer(poly_name='PEO')
	NPEO = 8
	peo.gen_pos(Nrepeats=NPEO, topology='linear')

	# generate PVA backbone
	pva=Polymer(poly_name='PVA')
	NPVA = 50
	topology='linear'
	pva.gen_pos(Nrepeats=NPVA, topology=topology)

	# generate the bottle brush PVA-g-PEO
	pva_g_peo = polyGraft(pva, peo)

	# set grafting density: grafts/monomer
	gft = 0.2
	pva_g_peo.setGraftingDensity(gft)

	# generate the graft structure
	pva_g_peo.setGftAtoms('O')
	pva_g_peo.genGraftStruct()

	# # save pdb and rtp files
	pva_g_peo.toPDB(f"{topology}-PVA{NPVA}_PEO{NPEO}_sigma_{gft}.pdb")
	pva_g_peo.toRTP(f"{topology}-PVA{NPVA}_PEO{NPEO}_sigma_{gft}.rtp")
```