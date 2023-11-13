This is the run directory to generate the structure and topology of polymer-grafted hybrid materials, with examples given below.
- poly-g-hard
  - nanoslab brush
  - nanoparticle brush
  - nanorod brush
  - nanopore brush

    
- poly-g-soft
  - bottlebrush polymer with a linear backbone
  - bottlebrush polymer with a cyclic backbone

Note: there are two versions of the generation code. Generation python file and folder name ending with *_lmp* (e.g., genNP_lmp.py and nanoparticle_brush_lmp) refers to the generated polyGraft structure of lammps data file. Note also that the example data file given for PEO (PEO20.data) is a coil rather than a rod structure which is not usable for high grafting density cases. One needs to prepare a rod structure for that end. The mass is not included in the data file, one can add this information in the input file of running a simulation. A file named **polyGraft.prm** is used to specify the force field parameters for atoms/bonds/angles/dihedrals/impropers.

Note: there are two steps in the generation of poly-g-soft. First, use **genBBP.py** to generate the pdb and rtp file, then use _gmx pdb2gmx_ (Gromacs) to get the itp and gro file. For now, the grafting density ($\sigma$) of bottlebrush polymer cannot be arbitrary. One must ensure the backbone can be divided into equal lengths (_i.e._, N<sub>bb</sub> * $\sigma$ is an integer and is evenly divisible to the spacing distance).

Note: for cyclic bottlebrush polymers, there is a further step on the generated itp file by Gromacs. Add the bonding information produced in the ***_to_add_cyclic.top** file (in the file, only bond is given; if angles/dihedrals/pairs are also needed, add them as well).

To accommodate for different side chains and/or backbones (than the PVA-g-PEO example provided), one needs to update the residue/atom/bond definition in the **rtp_define.py** file as well as the **gen_pos** function in ../polyGraft/polymer.py. 
