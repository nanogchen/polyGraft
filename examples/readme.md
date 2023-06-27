This is the run directory to generate the structure and topology of polymer-grafted hybrid materials, with examples given below.
- poly-g-hard
  - nanoslab brush
  - nanoparticle brush
  - nanorod brush
  - nanopore brush

    
- poly-g-soft
  - bottlebrush polymer

Note: there are two steps in the generation of poly-g-soft. First, use **genBBP.py** to generate the pdb (without Hydrogens) and rtp file, add hydrogens in _e.g._, Avogadro, then use **rtp_define.py** to wrap up the updated pdb file. For now, the grafting density of bottlebrush polymer cannot be arbitrary. One must ensure it can be divided into equal lengths (_i.e._, Nbb*sigma is an integer).

To accommodate for different side chains and/or backbones (than the PVA-g-PEO example provided), one needs to update the residue/atom/bond definition in the **rtp_define.py** file. 
