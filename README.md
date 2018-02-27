zmattools
=========

Tools to generate and manipulate z-matrices

Description
===========
zmatools is a package for generating and manipulate the z-matrix of a small molecules

Currently it has two tools:

1. Generating the z-matrix from a structure file
2. Producing a trajectory when each dihedral is sampled in a specific range, assuming the rest of the atoms are fixed.

Installation
============
To install this tool, clone or download the git repository. Change to the downloaded directory and install the software with

```
python setup.py install
```

Examples
========

To find the available commands

```
zmat -h
```

To find the help of individual commands

```
zmat generate -h
```

To generate a z-matrix to be used with Gaussian

```
zmat generate model0.mol2 --gaussian
```

To obtain verbose information on the z-matrix generation

```
zmat generate model0.mol2 --verbosity 2
```

To produce a dihedral trajectory in the range -180 to 180 in 15 degrees increment

```
zmat traj model0.mol2
```

the trajectories will be in xyz-formats in ``model0_dihed_traj4.xyz``,
``model0_dihed_traj5.xyz``etc.

To produce a dihedral trajectory in the range -180 to 180 in 30 degrees increment

```
zmat traj model0.mol2 -d 30
```
