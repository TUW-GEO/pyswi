Introduction
============
Soil water index (SWI) is calculated from surface soil moisture (SSM) as a function of time needed for infiltration and can be used as a simple approximation of root-zone conditions.
SWI method was originally developed at TU Wien (Wagner et al. 1999; https://www.sciencedirect.com/science/article/pii/S003442579900036X) and later reformulated recursively by Albergel et al. (2008; https://hess.copernicus.org/articles/12/1323/2008/) to allow for production of an operational product in near-real time.
SWI is calculated for a given time length *T* using the formula:

.. math::

    SWI_T(t_n) = SWI_T(t_n _- _1) + gain_T(t_n)(SSM(t_n) - (SWI_T(t_n _- _1))

Where T-value is a temporal length ruling the infiltration; depending on the soil characteristics it translates into different soil depths.

Pyswi allows for calculation of SWI from SSM datasets using an exponential filter.


Features
=========
Calculation of SWI from SSM daily images as well as time series, examples enclosed.