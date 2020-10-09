Introduction
============
Soil water index (SWI) is calculated from surface soil moisture (SSM) as a function of time needed for infiltration by means of an exponential filter and can be used as a simple approximation of root-zone conditions.
SWI method was originally developed at TU Vienna (Wagner et al. 1999) and later reformulated recursively by Albergel et al. (2008) to allow for production of an operational product in near-real time.

Temporal length *T* rules the infiltration; depending on the soil characteristics it translates into different soil depths.
SWI is calculated for a given *T-value* using the formula (where :math:`_n` refers to the latest, and :math:`_{n-1}` to the previous SSM observation):

.. math::

    SWI(t_n) = SWI_T(t_{n-1}) + gain_T(t_n)(SSM(t_n) - SWI_T(t_{n-1}))

:math:`gain_T` at time :math:`(t_n)` is calculated as follows:

.. math::

    gain_T(t_n) = \frac{gain_T(t_{n-1})}{gain_T(t_{n-1})+e^-\frac{t_n - t_{n-1}}{T}}

Pyswi allows for calculation of SWI from SSM datasets, the process can be initialized at any time when SSM observations become available with the following values: :math:`SWI_T(t_n) = SSM(t_n), gain_T(t_n) = 1, t_{n-1} = t_n`

Features
========
    * Calculation of SWI from SSM images with regular (*by default*), and irregular (parameter *weighted*; Bauer-Marschalinger and Paulik, 2019) timestamps.

    * Calculation of SWI from SSM time series.

Examples enclosed.

References
==========
Albergel, C., Rüdiger, C., Pellarin, T., Calvet, J.-C., Fritz, N., Froissard, F., Suquia, D., Petitpa, A., Piguet, B., and Martin, E.: From near-surface to root-zone soil moisture using an exponential filter: an assessment of the method based on in-situ observations and model simulations, Hydrol. Earth Syst. Sci., 12, 1323–1337, https://doi.org/10.5194/hess-12-1323-2008, 2008.

Bauer-Marschalinger, B. and Paulik, C. (2019) Soil Water Index Algorithm Theoretical Basis Document. Copernicus Global Land Operations: Vegetation and Energy, CGLOPS-1. l1.20, issued 24.04.2019.

Wagner, W., Lemoine, G., and Rott, H.: A method for estimating soil moisture from ERS scatterometer and soil data, Remote
Sens. Environ., 70, 191–207, https://doi.org/10.1016/S0034-4257(99)00036-X, 1999.