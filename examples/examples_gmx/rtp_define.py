import sys
import math
import MDAnalysis as mda

# residue definition dictionary of atoms
res_rtp_dict = {
	"C3": "opls_135   -0.18    1", # C in -CH3 of backbone
	"HE": "opls_140    0.06    1", # H in -CH3 of backbone
	"C2": "opls_136   -0.12    1", # C in -CH2- of backbone
	"HP": "opls_140    0.06    1", # H in -CH2- of backbone
	"CT": "opls_183    0.18    1", # C in -C(H)-OP/OH
	"HT": "opls_185    0.06    1", # H in -C(H)-OP/OH
	"OP": "opls_180   -0.48    1", # O in -OCC of PEO
	"CS": "opls_182    0.12    1", # C in -OCC of PEO
	"H" : "opls_185    0.06    1", # H in -OCC of PEO
	"OH": "opls_154   -0.68    1", # O in -OH
	"HO": "opls_155    0.44    1", # H in -OH
	"OA": "opls_154   -0.68    1", # O in -OH
	"HA": "opls_155    0.44    1"  # H in -OH
}

def gen_sc_rtp(polyname, Nrepeats):
	# generate the rtp files to be used for side-chain grafting

	BBP_sc_atoms = []
	BBP_sc_bonds = []

	if polyname == "PEO":
		for i in range(Nrepeats):
			BBP_sc_atoms.append(f'OP{i+1}')
			BBP_sc_atoms.append(f'CS{2*i+1}')
			BBP_sc_atoms.append(f'H{4*i+1}')
			BBP_sc_atoms.append(f'H{4*i+2}')
			BBP_sc_atoms.append(f'CS{2*i+2}')
			BBP_sc_atoms.append(f'H{4*i+3}')
			BBP_sc_atoms.append(f'H{4*i+4}')

			# bonds
			BBP_sc_bonds.append(f"OP{i+1} CS{2*i+1}")
			BBP_sc_bonds.append(f"CS{2*i+1} H{4*i+1}")
			BBP_sc_bonds.append(f"CS{2*i+1} H{4*i+2}")
			BBP_sc_bonds.append(f"CS{2*i+1} CS{2*i+2}")
			BBP_sc_bonds.append(f"CS{2*i+2} H{4*i+3}")
			BBP_sc_bonds.append(f"CS{2*i+2} H{4*i+4}")
			
			if i != Nrepeats-1:
				BBP_sc_bonds.append(f"CS{2*i+2} OP{i+2}")

		# add the tail
		BBP_sc_atoms.append("OA")
		BBP_sc_atoms.append("HA")
		BBP_sc_bonds.append(f"CS{2*i+2} OA")
		BBP_sc_bonds.append(f"OA HA")

	elif polyname == "your-polymer":
		print(f"Please define your polymer here! Exiting...")
		sys.exit(0)

	else:
		print(f"Unknown polymer type ({polyname}) other than PEO! Exiting...")
		sys.exit(0)

	return BBP_sc_atoms, BBP_sc_bonds

def gen_BBP_rtp(Nsc, spacing):

	BBP_G_atoms=[]
	BBP_G_bonds=[]
	BBP_sc_atoms, BBP_sc_bonds = gen_sc_rtp('PEO', Nsc)

	# add bb (PVA) first and then sc (PEO): first with grafts; others no
	# CH2
	BBP_G_atoms.append(f"C2{1}")
	BBP_G_atoms.append(f"HP{1}")
	BBP_G_atoms.append(f"HP{2}")
	BBP_G_bonds.append(f"C2{1} HP{1}")
	BBP_G_bonds.append(f"C2{1} HP{2}")

	# CH
	BBP_G_atoms.append(f"CT{1}")
	BBP_G_atoms.append(f"HT{1}")
	BBP_G_bonds.append(f"C2{1} CT{1}")
	BBP_G_bonds.append(f"CT{1} HT{1}")
	BBP_G_bonds.append(f"CT{1} {BBP_sc_atoms[0]}")

	# side chain
	BBP_G_atoms.extend(BBP_sc_atoms)
	BBP_G_bonds.extend(BBP_sc_bonds)

	# the other repeat units without grafts
	i = 1
	for i in range(1,spacing):
		# CH2
		BBP_G_atoms.append(f"C2{i+1}")
		BBP_G_atoms.append(f"HP{2*i+1}")
		BBP_G_atoms.append(f"HP{2*i+2}")
		BBP_G_bonds.append(f"CT{i} C2{i+1}")
		BBP_G_bonds.append(f"C2{i+1} HP{2*i+1}")
		BBP_G_bonds.append(f"C2{i+1} HP{2*i+2}")
		BBP_G_bonds.append(f"C2{i+1} CT{i+1}")

		# CH
		BBP_G_atoms.append(f"CT{i+1}")
		BBP_G_atoms.append(f"HT{i+1}")
		BBP_G_bonds.append(f"CT{i+1} HT{i+1}")

		# OH
		BBP_G_atoms.append(f"OH{i}")
		BBP_G_atoms.append(f"HO{i}")
		BBP_G_bonds.append(f"CT{i+1} OH{i}")
		BBP_G_bonds.append(f"OH{i} HO{i}")

	# add bonds to the previous/next
	if spacing == 1:
		BBP_G_bonds.append(f"-CT{i} C2{1}")
		BBP_G_bonds.append(f"CT{i} +C2{1}")
	else:
		BBP_G_bonds.append(f"-CT{i+1} C2{1}")
		BBP_G_bonds.append(f"CT{i+1} +C2{1}")
	
	return BBP_G_atoms, BBP_G_bonds

def findGft(Nbb,gft_density):
	# find the index of the repeat unit that has a graft
	
	assert gft_density<=1.0, "Grafting density should not be larger than 1.0!"
	spacer = int(1.0/gft_density)

	gft_idx = []
	for i in range(Nbb):

		# add side chain if satisfy the criteria
		if i % spacer == 0:
			gft_idx.append(i)

	return gft_idx

