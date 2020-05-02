=========
Changelog
=========

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
