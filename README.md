zmattools
=========

Tools to generate and manipulate z-matrices

Description
===========
zmatools is a package for generating a z-matrix of a small molecules

Currently this is the only "tool", but plans are to add more functionality

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
