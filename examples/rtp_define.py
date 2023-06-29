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
	"HO": "opls_155    0.44    1"  # H in -OH
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
		BBP_sc_atoms.append("OH")
		BBP_sc_atoms.append("HO")
		BBP_sc_bonds.append(f"CS{2*i+2} OH")
		BBP_sc_bonds.append(f"OH HO")

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

def PDBwrap(fname, raw_uni, Nsc, spacing, Ngrafts):
	# wrap up pdb files of updated molecular file (after added Hydrogen in Avogadro)

	BTB_G_atoms, _ = gen_BBP_rtp(Nsc, spacing)

	# get res for 1/2/3
	BTB2_G = BTB_G_atoms.copy()
	BTB1_G = ['C31','HE1','HE2','HE3']
	BTB1_G.extend(BTB_G_atoms)
	BTB3_G = BTB_G_atoms.copy()
	BTB3_G.extend(['C31','HE1','HE2','HE3'])

	# read in the updated pdb file
	C_pos = raw_uni.select_atoms('name C').positions.tolist()
	O_pos = raw_uni.select_atoms('name O').positions.tolist()
	H_pos = raw_uni.select_atoms('name H').positions.tolist()

	# use a dict to save the xyz
	all_pos = {'O':O_pos,\
			   'C':C_pos,\
			   'H':H_pos}

	# residues: BTB1 (1) + BTB2 (N-2) + BTB3 (1)	
	with open(fname, 'w') as FO:

		# header
		FO.write("TITLE Write by polyGraft-pdbwrap \n")

		# atomidx
		atomid = 0

		# out variables
		rectype = 'ATOM'

		# loop over residues
		for i_Gft in range(Ngrafts):
			resid = i_Gft+1

			# BTB1
			if i_Gft == 0:

				for atom_idx in range(len(BTB1_G)):
					atomid += 1
					atomname =  BTB1_G[atom_idx]
					resname = 'BTB1'
					xyz = 	  all_pos[BTB1_G[atom_idx][0]].pop(0)
					occ = 	  1.00
					tempfac = 0.00
					element = BTB1_G[atom_idx][0]

					FO.write("{:<6s}{:>5d} {:>4s} {:>3s} {:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>12s}\n".format(\
							rectype, atomid, atomname, resname, resid, xyz[0], xyz[1], xyz[2],occ,tempfac,element))

			# BTB3
			elif i_Gft == Ngrafts-1:

				for atom_idx in range(len(BTB3_G)):
					atomid += 1
					atomname =  BTB3_G[atom_idx]
					resname = 'BTB3'
					xyz = 	  all_pos[BTB3_G[atom_idx][0]].pop(0)
					occ = 	  1.00
					tempfac = 0.00
					element = BTB3_G[atom_idx][0]

					FO.write("{:<6s}{:>5d} {:>4s} {:>3s} {:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>12s}\n".format(\
							rectype, atomid, atomname, resname, resid, xyz[0], xyz[1], xyz[2],occ,tempfac,element))

			# BTB2
			else:

				for atom_idx in range(len(BTB2_G)):
					atomid += 1
					atomname =  BTB2_G[atom_idx]
					resname = 'BTB2'
					xyz = 	  all_pos[BTB2_G[atom_idx][0]].pop(0)
					occ = 	  1.00
					tempfac = 0.00
					element = BTB2_G[atom_idx]

					FO.write("{:<6s}{:>5d} {:>4s} {:>3s} {:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>12s}\n".format(\
							rectype, atomid, atomname, resname, resid, xyz[0], xyz[1], xyz[2],occ,tempfac,element))

		# end
		FO.write("END\n")

# ------------------------------------------------------------------		

