# ProLIF-Coloring
## Overview

This is an C++ implementation of algorithms capable of generate discrete volume-application-descriptions of protein-ligand interactions.

This implementation is based on [ProLIF](https://github.com/chemosim-lab/ProLIF/) protein-ligand interaction definitions

This implementation uses [RDKit](https://www.rdkit.org/) software library in order to handle SMART patterns, geometry calculation and .pdb files import/export

This is part of a project-work for Advanced Computer Architectures course [@Polimi](https://www.polimi.it/)

### Authors
Project developed by:
- [Enrico Tirri](https://github.com/EnricoTirri)

### File Structure

* `include-base` -  headers files where:
    * `Mesh.hpp` - defines the data structure containing the discrete molecule mesh
    * `Discretizer.hpp` - defines method to transform rdkit molecule into mesh and vice versa
    * `Interaction.hpp` - defines the basic method interfaces and structure an interaction should have
    * `InteractionCollection.hpp` - defines the structure the collection of interaction to be applied should have
    * `DistanceInteraction.hpp` - defines the basic method interface and structure a distance-based interaction should have
    * `SingleAngleInteraction.hpp` - defines the basic method interface and structure a single-angle-based interaction should have
    * `RestrictedBasePIStackingInteraction.hpp` - defines the basic method interface and structure a base-pi-stacking-based interaction should have

* `include-extended` - header files of interaction classes extensions

* `src` - source file of headers multi-platform implementation

* `src-normal` - source file of headers cpu based, single threaded implementation

* `src-omp` - source file of headers cpu based, openmp multithreaded implementation

* `src-cuda` - source file of headers gpu based, cuda implementation

### How to build

In order to build the executable, from the root folder run the following commands:

```bash
$ mkdir build
$ cd build
$ cmake .. _FLAGS_
$ make
```
`_FLAGS_` are optional, they can be:
* `-D USECUDA=1/0'` - specify if to use gpu based, cuda implementation of interactions
* `-D CUDA_BLOCK_SIZE=_size_block_` - specify the size of cuda-thread-block to be used
* `-D USEOMP=1/0'` - specify if to use cpu base, openmp implementation of interactions
* `-D OMP_NUM_THREADS=_num_threads_` - specify the number of threads used by omp implementation

### How to run

```bash
$ ./ProLIF_Coloring input.pdb
```
where:
* `input.pdb` - is a path/filename to .pdb input molecule file

### Showcase

| ![Molecule](showcase/mol.gif)                  | ![DiscreteMolecule](showcase/dicr_mol.gif) |
|------------------------------------------------|--------------------------------------------|
| ![HydrofobicInteraction](showcase/hyd_int.gif) |                                            |