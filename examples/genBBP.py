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
	# peo.polyGRO_.atoms.write(f"PEO{NPEO}.pdb")

	# generate PVA backbone
	pva=Polymer(poly_name='PVA')
	NPVA = 50
	pva.gen_pos(Nrepeats=NPVA, topology='linear')
	# pva.polyGRO_.atoms.write(f"PVA{NPVA}.pdb")

	# generate the bottle brush PVA-g-PEO
	pva_g_peo = polyGraft(pva, peo)

	# set grafting density: grafts/monomer
	gft = 0.25
	pva_g_peo.setGraftingDensity(gft)

	# generate the graft structure
	pva_g_peo.setGftAtoms('O')
	pva_g_peo.genGraftStruct()

	# # save pdb and rtp files
	pva_g_peo.toPDB(f"PVA{NPVA}_PEO{NPEO}_sigma_{gft}.pdb")
	# pva_g_peo.toRTP(f"PVA{NPVA}_PEO{NPEO}_sigma_{gft}.rtp")
