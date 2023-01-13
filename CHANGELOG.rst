=========
Changelog
=========

Unreleased Changes in master
============================

-


Version 1.0.0 (2023-01-13)
==========================

- Includes improved error estimation scheme from `Pasik et al. (2023, submitted)`
  and `De Santis and Biondi (2018) <https://doi.org/10.29007/kvhb>`_.
- Public release on pypi (``pip install pyswi``)
- Inclusion of iterative storage module

Version 0.6.0 (TUW internal / no public release)
================================================

- Merged calc.py and iterative_swi.py into a single file
- Weights now optional in the original code (weighted functions removed)
- Unified variable names
- Updated the environment for py3, dropped py2 support

Version 0.5.0 (TUW internal / no public release)
================================================

- Python compatibility >= 3.6
- Update pyscaffold v3.2.3

Version 0.4.0 (TUW internal / no public release)
================================================

- Fixing performance bug in swi_ts by avoiding pandas.to_datetime function
- New simplified interface for swi_ts

Version 0.3.5 (TUW internal / no public release)
================================================

- Adds processing_start/end to empty_data

Version 0.3.4 (TUW internal / no public release)
================================================

- Weighted SWI does not calculate QFLAG by default (which should be done only daily). Needs proper implementation in the future.
- Weighted SWI rectifies overshooting due to weighting (--> 0% > SWI < 100%)

Version 0.3.3 (TUW internal / no public release)
================================================

- Correct initialisation of weighted SWI

Version 0.3.2 (TUW internal / no public release)
================================================

- Remove future dependency

Version 0.3.1 (TUW internal / no public release)
================================================
- Fix wrong character
