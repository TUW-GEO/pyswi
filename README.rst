=====
pyswi
=====


Python package for the swi calculation.


Description
===========

Python package for the swi calculation.


Installation
===========

Needed python packages see requirments.txt

swi_calc_routines.pyx has to be compiled using the following commands (Linux):

cython -a swi_calc_routines.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o swi_calc_routines.so swi_calc_routines.c

Note
====

This project has been set up using PyScaffold 2.5.7. For details and usage
information on PyScaffold see http://pyscaffold.readthedocs.org/.
