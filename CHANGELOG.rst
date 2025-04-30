=========
Changelog
=========

Unreleased Changes in master
============================
-

Version 1.2.0 (2025-04-30)
==========================
- Updated the swi_error_prop function so it works when creating 10 day RZSM updates for C3S. Before, the function didn't work properly when the given SSM input data contained less than 2 data points.

Version 1.1.0 (2024-06-06)
==========================

- SWI calculation in swi_error_prop() decoupled form input uncertainty data availability
- quality flags normalized within swi_error_prop() and returned as percentages
- restarting calculation with stored parameters (gain_in/out)
- tests for above changes and swi_error_prop()/calc_swi_ts() result equivalence

Version 1.0.0 (2023-01-13)
==========================

- Includes an uncertainty characterization scheme from `Pasik et al. (2023, in review) <https://doi.org/10.5194/egusphere-2023-47>`_.
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
