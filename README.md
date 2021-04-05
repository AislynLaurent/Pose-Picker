# Pose-Picker

[![DOI](https://zenodo.org/badge/306729907.svg)](https://zenodo.org/badge/latestdoi/306729907)

A program designed to select poses produced by molecular dynamics (MD) simulations. This project takes advantage of the [MDTraj library](https://github.com/mdtraj/mdtraj) to analyze a wide range of outputs from different MD software, and [SciKitLearn's](https://github.com/mdtraj/mdtraj) **[Agglomerative Clustering](https://scikit-learn.org/stable/modules/clustering.html#hierarchical-clustering)** packages for machine learning.

### Requires Python 3.7.3

## Scope

This program is still in development, and we're still working hard to improve it. Right now it works well for simple systems, in particular enzymes with a single active site or with multiple active sites that preform the same function. Work is still being done to evaluate preformance on more complex systems.

## Usage

Run this like a python script in a terminal window and follow the prompts. **You will need:**
* A folder containing trajectory files
  * Supported file formats include: `pdb`, `xtc`, `trr`, `dcd`, `binpos`, `netcdf`, `mdcrd`, `prmtop`, `gsd`... and more
  * Find more information about supported formats [here](http://mdtraj.org/)
* A PDB file which describes the system as it appears in the first frame
  * This is for labelling purposes - a raw PDB may not reflect preprocessing you have done
  * Innacurate PDB files will not throw an error unless there is a mismatch in the number of atoms
* A list of atom types of interest
  * A sample list is included here (based on AMBER / NAMD atom names) but you may supply your own
* The residue name for your subtrate OR residue of interest
  * You may get better results from residues from your substrate, depending on the system

## Issue Reporting and Updates

If you have feedback or if you encounter an error, please feel free to report it using GitHubs **[Issues Tab](https://github.com/AislynLaurent/Pose-Picker/issues)**. We check this page regularly and will do our best to attend to errors as quickly as we can.

# Thank You!
