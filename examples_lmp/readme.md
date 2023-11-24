This is the run directory to generate the structure and topology of polymer-grafted hybrid materials, with examples given below.
- poly-g-hard
  - nanoslab brush
  - nanoparticle brush
  - nanorod brush
  - nanopore brush

Note: this is the LAMMPS version of the generation code. Note also that the example data file given for PEO (PEO20.data) is a coil rather than a rod structure which is not usable for high grafting density cases. One needs to prepare a rod structure for that end. The mass is not included in the data file, one can add this information in the input file of running a simulation. A file named **polyGraft.prm** is used to specify the force field parameters for atoms/bonds/angles/dihedrals/impropers.
