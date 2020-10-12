*****
pyswi
*****

Python package allowing computation of the Soil Water Index from surface soil moisture observations by means of exponential filter.

Description
===========

The package includes the following features:

* SWI time series calculation from SSM time series
    * including two methods for SWI noise calculation (regular and recursive)
* Recursive SWI approach to calculate SWI for a single or a set of T-values in near-real time
    * also *Weighted* calculation of the SWI, allowing for custom weight assignment to individual observations

Installation
============
This package should be installable through pip:

.. code-block:: python

    pip install pyswi

Note
====

This project has been set up using PyScaffold 3.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
