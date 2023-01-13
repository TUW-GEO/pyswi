=====
pyswi
=====

Python package allowing computation of the Soil Water Index from surface soil moisture observations by means of exponential filter.

Description
===========

The package includes the following features:

* SWI time series calculation and error propagation from SSM time series
    * a recursive approach to SWI and its noise with calculation routines in Cython
    * an equivalent SWI calculation in Python with an exponential-filter-based
      error propagation scheme
* Recursive SWI approach to calculate SWI for a single or a set of T-values in near-real time
    * also *Weighted* calculation of the SWI, allowing for custom weight assignment to
      individual observations

iterative_storage
=================

Storage of iteration data between processing runs.

Description
===========

In a process that works iteratively and needs to store some data
between processing runs, the classes in this package can be used to store
that data as netCDF files of any format. The main functionality of this package
is in the building of the storage filenames and in reading the correct iteration
data from the disk when the process is started again.


Installation
============
This package should be installable through pip:

.. code-block:: python

    pip install pyswi

Note
====

This project has been set up using PyScaffold 3.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
