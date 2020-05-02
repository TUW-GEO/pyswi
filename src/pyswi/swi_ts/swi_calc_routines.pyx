# Copyright (c) 2019, TU Wien, Department of Geodesy and Geoinformation (GEO).
# All rights reserved.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL VIENNA UNIVERSITY OF TECHNOLOGY,
# DEPARTMENT OF GEODESY AND GEOINFORMATION BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


from libc.math cimport exp
import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def swi_calc_cy(np.ndarray[np.double_t] juldate,
                np.ndarray[np.float32_t] ssm,
                np.ndarray[np.int_t] ctime,
                np.ndarray[np.double_t] swi_jd,
                np.ndarray[np.float_t] nom,
                np.ndarray[np.float_t] denom,
                double last_jd_var,
                np.ndarray[np.float_t] norm_factor,
                double nan):

    cdef int i
    cdef int j = 0
    cdef int k
    cdef np.ndarray[np.float_t] ef = np.empty(len(ctime))
    cdef double tdiff
    cdef float qflag
    cdef int len_swi = len(swi_jd)
    cdef int len_jd = len(juldate)
    cdef int len_ctime = len(ctime)
    cdef np.ndarray[np.float32_t, ndim = 2] swi_ts = np.empty((len_swi, len_ctime), dtype=np.float32)
    cdef np.ndarray[np.float32_t, ndim = 2] swi_qflag = np.empty((len_swi, len_ctime), dtype=np.float32)

    for i in range(len_swi):

        while j < len_jd and juldate[j] <= swi_jd[i]:
            if ssm[j] != nan or ssm[j] == np.nan:
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
            if ssm[i] != nan or ssm[i] == np.nan:
                swi_ts[i, k] = nom[k] / denom[k]
            else:
                swi_ts[i, k] = np.nan

            qflag = 100.0 * \
                denom[k] * exp(-(swi_jd[i] - last_jd_var) /
                               ctime[k]) / norm_factor[k]
            if qflag > 100:
                qflag = 100
            swi_qflag[i, k] = qflag

    return swi_ts, swi_qflag, nom, denom, last_jd_var


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def swi_calc_cy_noise(np.ndarray[np.double_t] juldate,
                      np.ndarray[np.float32_t] ssm,
                      np.ndarray[np.int_t] ctime,
                      np.ndarray[np.double_t] swi_jd,
                      np.ndarray[np.float_t] nom,
                      np.ndarray[np.float_t] denom,
                      double last_jd_var,
                      np.ndarray[np.float32_t] ssm_noise,
                      np.ndarray[np.double_t] nom_noise):

    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    cdef np.ndarray[np.float_t] ef = np.empty(len(ctime))
    cdef double tdiff = 0.0
    cdef int len_swi = len(swi_jd)
    cdef int len_jd = len(juldate)
    cdef int len_ctime = len(ctime)
    cdef np.ndarray[np.float32_t, ndim = 2] swi_ts = np.empty((len_swi, len_ctime), dtype=np.float32)
    cdef np.ndarray[np.float32_t, ndim = 2] swi_noise = np.empty((len_swi, len_ctime), dtype=np.float32)

    for i in range(len_swi):

        while j < len_jd and juldate[j] <= swi_jd[i]:
            tdiff = juldate[j] - last_jd_var
            for k in range(len_ctime):
                ef[k] = exp(-tdiff / ctime[k])
                nom[k] = ef[k] * nom[k] + ssm[j]
                denom[k] = ef[k] * denom[k] + 1
                exp_term = ef[k] ** 2
                nom_noise[k] = nom_noise[k] * exp_term + ssm_noise[j]
            last_jd_var = juldate[j]
            j += 1

        if denom[0] == 0:
            continue  # no valid SSM measurement before swi_jd[i]

        for k in range(len_ctime):
            swi_ts[i, k] = nom[k] / denom[k]
            denom_ns = denom[k] ** 2
            swi_noise[i, k] = nom_noise[k] / denom_ns

    return swi_ts, swi_noise, nom, denom, last_jd_var, nom_noise
