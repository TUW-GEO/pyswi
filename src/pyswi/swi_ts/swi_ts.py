# Copyright (c) 2021, TU Wien, Department of Geodesy and Geoinformation (GEO).
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

from math import exp

import pyximport
import numpy as np

pyximport.install(setup_args={'include_dirs': [
                  np.get_include()]}, inplace=True)

from pyswi.swi_ts.swi_calc_routines import swi_calc_cy, swi_calc_cy_noise

def calc_swi_ts(ssm_ts, swi_jd, gain_in=None, t_value=[1, 5, 10, 15, 20],
                nom_init=0, denom_init=0, nan=-9999.):
    """
    Time series calculation of the Soil Water Index.

    Parameters
    ----------
    ssm_ts : numpy.ndarray or dict
        Surface soil moisture time series with fields: jd, sm, sm_noise
    swi_jd : numpy.ndarray
        Julian date time stamps of the SWI time series.
    gain_in : dict, optional
        Gain parameters of last calculation.
        Dictionary with fields: last_jd, nom, denom
    t_value : list, optional
        Characteristic time length (default: 1, 5, 10, 15, 20).
    nom_init : float64, optional
        Initial value of nom in the SWI calculation (default: 0).
    denom_init : float64, optional
        Initial value of denom in the SWI calculation (default: 0).
    nan : float64, optional
        NaN value to be masked in the SWI retrieval.

    Returns
    -------
    swi_ts : numpy.ndarray
        Soil Water Index (SWI) time series.
    gain_out : dict
        Gain parameters of last calculation.
        fields gpi, last_jd, nom, denom, nom_ns
    """

    t_value = np.asarray(t_value)

    if gain_in is None:
        gain_in = {'last_jd': np.float64(0),
                   'denom': np.full(len(t_value), denom_init, dtype=np.float64),
                   'nom': np.full(len(t_value), nom_init, dtype=np.float64),
                   'nom_noise': np.zeros(len(t_value))}

    if gain_in['denom'][0] == 1:
        last_jd_var = ssm_ts['jd'][0]
    else:
        last_jd_var = gain_in['last_jd']

    norm_factor = \
        np.sum(np.exp(-(np.tensordot(np.arange(1000, dtype=np.float32),
                                     (1. / t_value), axes=0))), axis=0)

    nom_noise = np.ones(len(t_value), dtype=np.float64)

    if 'sm_noise' in ssm_ts.dtype.fields:
        if gain_in['nom_noise'] is not None:
            nom_noise = gain_in['nom_noise']

        swi, swi_noise, nom, denom, last_jd_var, nom_ns = \
            swi_calc_cy_noise(ssm_ts['jd'], ssm_ts['sm'], t_value, swi_jd,
                              gain_in['nom'], gain_in['denom'], last_jd_var,
                              ssm_ts['sm_noise'], nom_noise)
        swi_qflag = np.zeros_like(swi)
    else:
        swi, swi_qflag, nom, denom, last_jd_var = \
            swi_calc_cy(ssm_ts['jd'], ssm_ts['sm'], t_value, swi_jd,
                        gain_in['nom'], gain_in['denom'], last_jd_var,
                        norm_factor, nan)
        swi_noise = np.zeros_like(swi)

    gain_out = {'denom': denom, 'nom': nom,
                'last_jd': last_jd_var, 'nom_noise': nom_noise}

    dtype_list = [('jd', np.float64)]

    for t in t_value:
        dtype_list.append(('swi_{}'.format(t), np.float32))
        dtype_list.append(('swi_noise_{}'.format(t), np.float32))
        dtype_list.append(('swi_qflag_{}'.format(t), np.float32))

    swi_ts = np.zeros(swi_jd.size, dtype=np.dtype(dtype_list))

    swi_ts['jd'] = swi_jd
    for i, t in enumerate(t_value):
        swi_ts['swi_{}'.format(t)] = swi[:, i]
        swi_ts['swi_noise_{}'.format(t)] = swi_noise[:, i]
        swi_ts['swi_qflag_{}'.format(t)] = swi_qflag[:, i]

    return swi_ts, gain_out


def calc_swi_noise(ssm_ts, t_value):
    """
    Calculation of Soil Water Index (SWI) noise.
    (Slower than the recursive approach).

    Parameters
    ----------
    ssm_ts : numpy.ndarray
        Surface soil moisture time series with fields: jd, sm, sm_noise
    t_value : numpy.ndarray
        Characteristic time length.

    Returns
    -------
    swi_noise_ts : numpy.ndarray
        Soil Water Index noise time series.
    """
    len_sm = len(ssm_ts)
    len_t_value = len(t_value)

    nom_db = np.zeros((len_sm, len_t_value))
    den_db = np.zeros((len_sm, len_t_value))
    swi_noise = np.zeros((len_sm, len_t_value))

    for i in range(0, len_sm):
        nom = np.zeros(len_t_value)
        den = np.zeros(len_t_value)
        for c in range(0, len_t_value):
            for j in range(0, i+1):
                jd_diff = ssm_ts['jd'][i] - ssm_ts['jd'][j]
                exp_term = exp(((-2)*(jd_diff)) / t_value[c])
                nom[c] += ssm_ts['sm_noise'][i] * exp_term
                den[c] += exp(((-1)*(jd_diff)) / t_value[c])

            nom_db[i][c] = nom[c]
            den_db[i][c] = den[c]
            swi_noise[i][c] = (nom[c] / (den[c] ** 2))

    gain_out = {'denom': den, 'nom': nom, 'last_jd': ssm_ts['jd'][-1],
                'nom_noise': 0}

    dtype_list = [('jd', np.float64)]

    for t in t_value:
        dtype_list.append(('swi_noise_{}'.format(t), np.float32))

    swi_ts = np.zeros(ssm_ts.size, dtype=np.dtype(dtype_list))

    swi_ts['jd'] = ssm_ts['jd']

    for i, t in enumerate(t_value):
        swi_ts['swi_noise_{}'.format(t)] = swi_noise[:, i]

    return swi_ts, gain_out


def calc_swi_noise_rec(ssm_ts, t_value, last_den=1, last_nom=0):
    """
    Recursive calculation of Soil Water Index (SWI) noise.

    Parameters
    ----------
    ssm_ts : numpy.ndarray
        Surface soil moisture time series with fields: jd, sm, sm_noise
    t_value : numpy.ndarray
        Characteristic time length.
    denom : float
        denom value of the last calculation and starting point for
        the calculation.
    nom : float
        nom value of the last calculation and starting point for
        the calculation.

    Returns
    -------
    swi_noise_ts : numpy.ndarray
        Soil Water Index noise time series.
    """
    len_sm = len(ssm_ts)
    len_t_value = len(t_value)

    # var{n_swi} calculation
    nom_sn = np.zeros((len_sm, len_t_value))
    den = np.zeros((len_sm, len_t_value))
    swi_noise = np.zeros((len_sm, len_t_value))

    last_jd = ssm_ts['jd'][0]
    den_values = np.ones(len_t_value)
    nom_values = np.zeros(len_t_value)

    den_values.fill(last_den)
    nom_values.fill(last_nom)

    for i in range(0, len_sm):
        for c in range(0, len_t_value):
            tdiff = (ssm_ts['jd'][i]-last_jd) / t_value[c]
            fn = exp((-1) * tdiff)
            den[i][c] = 1 + fn * den_values[c]
            den_sn = den[i][c] ** 2

            exp_term = exp((-2)*tdiff)
            nom_sn[i][c] = nom_values[c] * exp_term + ssm_ts['sm_noise'][i]

            swi_noise[i][c] = nom_sn[i][c] / den_sn
            den_values[c] = den[i][c]
            nom_values[c] = nom_sn[i][c]
        last_jd = ssm_ts['jd'][i]

    gain_out = {'denom': den_values, 'nom': nom_values, 'last_jd': last_jd,
                'nom_noise': 0}

    dtype_list = [('jd', np.float64)]
    for t in t_value:
        dtype_list.append(('swi_noise_{}'.format(t), np.float32))

    swi_ts = np.zeros(ssm_ts.size, dtype=np.dtype(dtype_list))

    swi_ts['jd'] = ssm_ts['jd']
    for i, t in enumerate(t_value):
        swi_ts['swi_noise_{}'.format(t)] = swi_noise[:, i]

    return swi_ts, gain_out