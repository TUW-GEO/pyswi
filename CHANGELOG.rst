=========
Changelog
=========

Version 0.6.0
=============
- Merged calc.py and iterative_swi.py into a single file
- Weights now optional in the original code (weighted functions removed)
- Unified variable names
- Updated the environment for py3, dropped py2 support

Version 0.5.0
=============

- Python compatibility >= 3.6
- Update pyscaffold v3.2.3

Version 0.4.0
=============

- Fixing performance bug in swi_ts by avoiding pandas.to_datetime function
- New simplified interface for swi_ts

Version 0.3.5
=============

- Adds processing_start/end to empty_data

Version 0.3.4
=============

- Weighted SWI does not calculate QFLAG by default (which should be done only daily). Needs proper implementation in the future.
- Weighted SWI rectifies overshooting due to weighting (--> 0% > SWI < 100%)

Version 0.3.3
=============

- Correct initialisation of weighted SWI

Version 0.3.2
=============

- Remove future dependency

Version 0.3.1
=============

- Fix wrong character
