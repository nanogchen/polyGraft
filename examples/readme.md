This is the run directory to generate the structure and topology of polymer-grafted hybrid materials, with examples given below.
- poly-g-hard
  - nanoslab brush
  - nanoparticle brush
  - nanorod brush
  - nanopore brush

    
- poly-g-soft
  - bottlebrush polymer with linear backbone
  - bottlebrush polymer with cyclic backbone


Note: there are two steps in the generation of poly-g-soft. First, use **genBBP.py** to generate the pdb (without Hydrogens) and rtp file, add hydrogens in _e.g._, Avogadro (or use [openbabel](https://openbabel.org/docs/dev/Command-line_tools/babel.html), see command below), then use **rtp_define.py** to wrap up the updated pdb file, such that one can use _gmx pdb2gmx_ to get the itp and gro file. For now, the grafting density of bottlebrush polymer cannot be arbitrary. One must ensure it can be divided into equal lengths (_i.e._, Nbb*sigma is an integer).

```
obabel -ipdb onlyHeavyAtoms.pdb -O allAtoms_wHs.pdb -h
```

To accommodate for different side chains and/or backbones (than the PVA-g-PEO example provided), one needs to update the residue/atom/bond definition in the **rtp_define.py** file as well as the **gen_pos** function in ../polyGraft/polymer.py. 
