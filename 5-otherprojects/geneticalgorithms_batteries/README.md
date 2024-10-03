# Random Generation and Genetic Algorithm for Bipyridine-Type Molecules

This project contains two components that work together to randomly generate 
bipyridine (bpy)-type molecules and then apply a genetic algorithm (GA) to 
optimize the generated molecules based on functional groups, double bonds, 
and nitrogen positions. The two core scripts are `run-random-gen.ipynb` and 
`run-GA.ipynb`, each serving a distinct purpose in the workflow.

---

## Random Generation of Bipyridine-Type Molecules

### `run-random-gen.ipynb`

This script performs the random generation of bpy-type molecules by adding 
various elements such as functional groups (FGs), double bonds (DBs), and 
nitrogen atoms (N) in different positions around the bipyridine ring. The 
script generates combinatorial possibilities by randomly selecting molecular 
features and placing them symmetrically around the ring.

#### Key Arguments:

- **`INITIAL_MOL_FILEPATH`**: Path to the starting molecule file (in `.mol` 
  format) that serves as the base structure for generating new molecules.
- **`MAX_NUMBER_N`**: Maximum number of nitrogen atoms allowed to be added 
  to the ring.
- **`MAX_NUMBER_FG`**: Maximum number of functional groups allowed to be 
  added to the ring.
- **`MAX_NUMBER_DB`**: Maximum number of double bonds allowed to be added 
  to the ring.
- **`MOLECULE_NAME`**: Name of each generated molecule, usually formatted 
  with an iteration number.
- **`NUMBER_TO_GENERATE`**: The number of molecules to randomly generate.

#### Example of Inputs:

```python
INITIAL_MOL_FILEPATH = './buildingblocks/startingmolfile.mol'
MAX_NUMBER_N = 4
MAX_NUMBER_FG = 2
MAX_NUMBER_DB = 6
DOUBLE = True  # Enable a double linker between rings
NUMBER_TO_GENERATE = 100  # Generate 100 molecules

# Genetic Algorithm for Molecule Optimization

### `run-random-gen.ipynb`

This project implements a genetic algorithm to optimize bipyridine-type 
(bpy-type) molecules. The algorithm combines two parent molecules and 
generates child molecules, potentially improving their properties based 
on genetic crossover operations. The goal is to fine-tune molecular 
structures by adding nitrogen atoms, functional groups, and double bonds 
to specific positions on the molecule.

---

## Overview

The genetic algorithm uses a starting bpy molecule and allows for various 
genetic modifications (crossover, mutation) to produce child molecules with 
desired properties. The core functionality includes performing crossover 
operations between two molecules, symmetrizing the structures, adding 
functional groups and bonds, and refining the final structure.

---

## Key Arguments

- **`INITIAL_MOL_FILEPATH`**: Path to the initial bpy molecule file that serves 
  as the starting point for genetic crossover.
- **`MAX_NUMBER_N`**: Maximum number of nitrogen atoms that can be added after 
  the genetic crossover operation.
- **`MAX_NUMBER_FG`**: Maximum number of functional groups to add after 
  crossover.
- **`MAX_NUMBER_DB`**: Maximum number of double bonds to add after crossover.
- **`DOUBLE`**: A boolean flag indicating whether a double bond linker is 
  involved in the genetic crossover process.

---

## Example of Inputs

Here’s how you can define your inputs for running the genetic algorithm:

```python
INITIAL_MOL_FILEPATH = './buildingblocks_GA/startingmolfile.mol'
MAX_NUMBER_N = 4
MAX_NUMBER_FG = 2
MAX_NUMBER_DB = 3
DOUBLE = True  # Enable a double linker during crossover

