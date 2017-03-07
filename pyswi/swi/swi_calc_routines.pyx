'''
Created on Oct 14, 2013

@author: Christoph Paulik chrisotph.paulik@geo.tuwien.ac.at
'''


cimport numpy as np
import numpy as np
from libc.math cimport exp
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def swi_calc_cy(np.ndarray[np.double_t] juldate,
                np.ndarray[np.float32_t] ssm,
                np.ndarray[np.int_t] ctime,
                np.ndarray[np.double_t] swi_jd,
                np.ndarray[np.float32_t, ndim=2] swi_ts,
                np.ndarray[np.float32_t, ndim=2] swi_qflag,
                np.ndarray[np.float_t] nom,
                np.ndarray[np.float_t] denom,
                double last_jd_var,
                np.ndarray[np.float_t] norm_factor):

    cdef int i
    cdef int j = 0
    cdef int k
    cdef np.ndarray[np.float_t] ef = np.empty(len(ctime))
    cdef double tdiff
    cdef float qflag
    cdef int len_swi = len(swi_jd)
    cdef int len_jd = len(juldate)
    cdef int len_ctime = len(ctime)

    for i in range(len_swi):

            while j < len_jd and juldate[j] <= swi_jd[i]:
                tdiff = juldate[j] - last_jd_var
                for k in range(len_ctime):
                    ef[k] = exp(-tdiff / ctime[k])
                    nom[k] = ef[k] * nom[k] + ssm[j]
                    denom[k] = ef[k] * denom[k] + 1

                last_jd_var = juldate[j]
                j += 1

            if denom[0] == 0:
                continue  # no valid SSM measurement before swi_jd[i]

            for k in range(len_ctime):

                swi_ts[i, k] = nom[k] / denom[k]
                qflag = 100.0 * denom[k] * exp(-(swi_jd[i] - last_jd_var) / ctime[k]) / norm_factor[k]
                if qflag > 100:
                    qflag = 100
                swi_qflag[i, k] = qflag

    return swi_ts, swi_qflag, nom, denom, last_jd_var
