This is the run directory to generate the structure and topology of polymer-grafted hybrid materials, with examples given below.
- poly-g-hard
  - nanoslab brush
  - nanoparticle brush
  - nanorod brush
  - nanopore brush

    
- poly-g-soft
  - bottlebrush polymer with linear backbone
  - bottlebrush polymer with cyclic backbone


Note: there are two steps in the generation of poly-g-soft. First, use **genBBP.py** to generate the pdb and rtp file, then use _gmx pdb2gmx_ (Gromacs) to get the itp and gro file. For now, the grafting density ($\sigma$) of bottlebrush polymer cannot be arbitrary. One must ensure the backbone can be divided into equal lengths (_i.e._, N<sub>bb</sub> * $\sigma$ is an integer).

Note: for cyclic bottlebrush polymers, there is a further step on the generated itp file by Gromacs. Add the bonding information produced in the ***_to_add_cyclic.top** file (in the file, only bond is given; if angles/dihedrals/pairs are also needed, add them as well).

To accommodate for different side chains and/or backbones (than the PVA-g-PEO example provided), one needs to update the residue/atom/bond definition in the **rtp_define.py** file as well as the **gen_pos** function in ../polyGraft/polymer.py. 
