=====
pyswi
=====

Python package for the swi calculation for time series and for images.

Description
===========

The package includes the following features:

* SWI TS: SWI TS Calculation with SWI GAIN support.
* SWI TS: Calculation of the QFLAG.
* SWI TS: Calculation of SWI NOISE including SWI NOISE GAIN.
* SWI IMG: Iterative approach to calculate SWI for a single T value.
* SWI IMG: Iterative approach to calculate SWI for a set of T values.

Installation
===========

In the root directory of the repository run the following command to create a
new conda environment:

conda env create -f conda_environment.yml

The new environment can be activated by calling:

source activate swi_env

Usage
===========

For the usage of the functions, please have a look at the unit tests located
in the tests folder.
