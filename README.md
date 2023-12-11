# polyGraft: A Program for Molecular Structure and Topology Generation of Polymer-Grafted Hybrid Nanostructures

- Unigrafts grafting in atomistic format (GROMACS and LAMMPS)
- Binary grafts in atomistic format (GROMACS and LAMMPS)
- Unigrafts grafting in coarse-grained format (LAMMPS)
- Binary grafts in coarse-grained format (LAMMPS)

<img src="polyGraft.jpg" alt="drawing" width="600"/>

# Prerequisites and installation
- [atomsk](https://atomsk.univ-lille.fr/) for hard substrate generation (any other tools should also be fine, save in pdb file format) (tested v0.12)
- [MDAnalysis](https://www.mdanalysis.org/) file IO (tested v2.5.0)
- [numba](https://numba.pydata.org/) accelerated array processing (tested v0.57.1)

Installation steps (tested with Anaconda):
1. install [Anaconda](https://anaconda.org/)
2. create a new environment adding MDAnalysis, numba etc, using the following command under the default env (base):
> conda env create -f environment.yml
3. download atomsk and place under /path/to/polyGraft/bin
4. download polyGraft through github
5. activate the polyGraft environment and use polyGraft for generation under /path/to/polyGraft/examples
> conda activate polyGraft

**NOTE**: to use other polymers (rather than PEO) for polymer brush generation using polyGraft, one must provide the gro/itp files (for Gromacs simulations) or data file (for LAMMPS simulations) of the polymer. The user should be responsible for the results in this case as no validity checks on the input files were applied (i.e., validate the connectivity or topology, and force field parameters used for further MD simulation etc.)

# How to cite
1. Chen, Guang. "polyGraft 1.0: A Program for Molecular Structure and Topology Generation of Polymer-Grafted Hybrid Nanostructures". J. Comput. Chem. 2023, 44(28), 2230. https://doi.org/10.1002/jcc.27206
2. Chen, Guang, and Elena E. Dormidontova. "Cyclic vs Linear Bottlebrush Polymers in Solution: Side-Chain Length Effect." Macromolecules 56.9 (2023): 3286–3295. https://doi.org/10.1021/acs.macromol.3c00362
3. Chen, Guang, and Elena Dormidontova. "PEO-Grafted Gold Nanopore: Grafting Density, Chain Length, and Curvature Effects." Macromolecules 55.12 (2022): 5222-5232. https://doi.org/10.1021/acs.macromol.2c00323

# How to use
Representative examples are given in the **examples*** directories, specifically
- *lmp: for lammps format
- *bi: for binary grafts used for grafting procedure
- *cg: coarse-grained data (bead-spring model)

For any nanostructure generation, it takes two steps generally:
1. generate or import the substrate material (hard or soft), and polymer (gro and itp or data files);
2. generate the structure and topology using polyGraft. For poly-g-hard, the files are in .gro/.itp or .data format, while for poly-g-soft, they are in pdb/rtp format.

# Documentation and User Guide
Please refer to the documentation of the code [wiki](https://github.com/nanogchen/polyGraft/wiki)!

# Seek help or new features
Open an [issue](https://github.com/nanogchen/polyGraft/issues)!
