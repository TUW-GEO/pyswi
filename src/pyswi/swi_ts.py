# Copyright (c) 2025, TU Wien.
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

from math import exp, sqrt

import pyximport
import numpy as np

pyximport.install(setup_args={'include_dirs': [np.get_include()]},
                  inplace=True)

from pyswi.swi_calc_routines import swi_calc_cy, swi_calc_cy_noise


def swi_error_prop(ssm, t_value, t_noise, swi_error, gain_in=None, nan=-9999.):
    """
    Recursive SWI calculation and error propagation function
    based on DeSantis and Biondi (2018; https://doi.org/10.29007/kvhb)
    Translated from MatLab code obtained from the authors.

    Parameters
    ----------
    ssm: numpy.ndarray
        Surface soil moisture time series with fields 'sm', 'sm_uncertainty' and 'sm_jd'
    t_value: numpy.ndarray
        Exponential filter characteristic T-value parameter
    t_noise: numpy.ndarray
        T-value standard error.
        10% of T for calibrated T-values.
    swi_error: numpy.ndarray
        Exponential filter model structural error.
        ubMSE(ISMNswi, ISMNrzsm), based on empirical experiments.
    gain_in: dict
        stored parameters of the last iteration.
    nan: float
        nan value of the input ssm dataset

    Returns
    -------
    swi: numpy.ndarray
        Soil water index time series
    swi_noise:numpy.ndarray
        Soil water index noise time series
    """

    len_ssm = len(ssm)
    len_T = len(t_value)
    swi = np.full((len_ssm, len_T), np.nan)
    swi_noise = np.full((len_ssm, len_T), np.nan)
    qflag = np.full((len_ssm, len_T), np.nan)
    qflag_norm = [
        np.sum(np.exp(-(np.tensordot(
            np.arange(1000, dtype=np.float32), (1. / t), axes=0))),
               axis=0) for t in t_value
    ]

    first_ssm_idx = np.argmax((ssm['sm'] != nan) & ~np.isnan(ssm['sm']))
    first_noise_idx = np.argmax((ssm['sm_uncertainty'] != nan)
                                & ~np.isnan(ssm['sm_uncertainty']))
    j = 1 + first_ssm_idx

    if gain_in is None:
        swi[first_ssm_idx] = ssm['sm'][first_ssm_idx]
        last_jd = ssm['sm_jd'][first_ssm_idx]
        gain_curr = [1] * len_T
        contr1_curr = [ssm['sm_uncertainty'][first_noise_idx]**2
                       ] * len_T  # delta
        G_curr = [0] * len_T
        JT_curr = [0] * len_T  # Jacobian term
        qflag[first_ssm_idx] = [1] * len_T
    else:
        time_diff = ssm['sm_jd'][first_ssm_idx] - gain_in['last_jd']
        ef = [np.exp(-time_diff / t_value[i]) for i in range(0, len(t_value))]
        gain_curr = [
            gain_in['last_gain'][i] / (gain_in['last_gain'][i] + ef[i])
            for i in range(0, len(t_value))
        ]
        swi[first_ssm_idx] = gain_in['last_swi'] + gain_curr * (
            ssm['sm'][first_ssm_idx] - gain_in['last_swi'])
        contr1_curr = [
            ((1 - gain_curr[i])**2) * gain_in['last_contr1'][i] +
            (gain_curr[i] * ssm['sm_uncertainty'][first_noise_idx])**2
            for i in range(0, len(t_value))
        ]
        G_curr = [
            ef[i] * (gain_in['last_G'][i] +
                     time_diff / t_value[i] / gain_in['last_gain'][i])
            for i in range(0, len(t_value))
        ]
        JT_curr = [
            gain_curr[i] / t_value[i] *
            (G_curr[i] *
             (gain_in['last_swi'][i] - swi[first_ssm_idx][i]) + ef[i] *
             t_value[i] / gain_in['last_gain'][i] * gain_in['last_JT'][i])
            for i in range(0, len(t_value))
        ]
        contr2 = [(JT_curr[i] * t_noise[i])**2 for i in range(0, len(t_value))]
        swi_noise[first_noise_idx] = [
            sqrt(contr1_curr[i] + contr2[i] + swi_error[i])
            for i in range(0, len(t_value))
        ]
        if ssm['sm'][0] != nan and ~np.isnan(ssm['sm'][0]):
            qflag[0] = [
                1 +
                gain_in['last_qflag'][i] * np.exp(-(1. / float(t_value[i])))
                for i in range(0, len(t_value))
            ]
        else:
            qflag[0] = [
                gain_in['last_qflag'][i] * np.exp(-(1. / float(t_value[i])))
                for i in range(0, len(t_value))
            ]
        for f in range(1, j):
            if ssm['sm'][f] != nan and ~np.isnan(ssm['sm'][f]):
                qflag[f] = [
                    1 + qflag[f - 1][i] * np.exp(-(1. / float(t_value[i])))
                    for i in range(0, len(t_value))
                ]
            else:
                qflag[f] = [
                    qflag[f - 1][i] * np.exp(-(1. / float(t_value[i])))
                    for i in range(0, len(t_value))
                ]
        last_jd = ssm['sm_jd'][first_ssm_idx]

    gain_old = [None] * len_T
    contr1_old = [None] * len_T
    contr2 = [None] * len_T
    G_old = [None] * len_T
    JT_old = [None] * len_T

    for i in range(j, len_ssm):
        while j < len_ssm and ssm['sm_jd'][j] <= ssm['sm_jd'][i]:
            if ssm['sm'][j] != nan and ~np.isnan(ssm['sm'][j]):
                time_diff = ssm['sm_jd'][j] - last_jd
                for c in range(len_T):
                    ef = np.exp(-time_diff / t_value[c])
                    gain_old[c] = gain_curr[c]
                    gain_curr[c] = gain_old[c] / (gain_old[c] + ef)
                    swi[i][c] = swi[i - int(time_diff)][c] + gain_curr[c] * (
                        ssm['sm'][i] - swi[i - int(time_diff)][c])
                    if ssm['sm_uncertainty'][j] != nan and ~np.isnan(
                            ssm['sm_uncertainty'][j]):
                        contr1_old[c] = contr1_curr[c]
                        contr1_curr[c] = (
                            (1 - gain_curr[c])**2) * contr1_old[c] + (
                                gain_curr[c] * ssm['sm_uncertainty'][i])**2
                        G_old[c] = G_curr[c]
                        JT_old[c] = JT_curr[c]
                        G_curr[c] = ef * (G_old[c] +
                                          time_diff / t_value[c] / gain_old[c])
                        JT_curr[c] = gain_curr[c] / t_value[c] * (
                            G_curr[c] *
                            (swi[i - int(time_diff)][c] - swi[i][c]) +
                            ef * t_value[c] / gain_old[c] * JT_old[c])
                        contr2[c] = (JT_curr[c] * t_noise[c])**2
                        swi_noise[i][c] = sqrt(contr1_curr[c] + contr2[c] +
                                               swi_error[c])
                last_jd = ssm['sm_jd'][i]
                last_valid_index = i
            j += 1

        for c in range(len_T):
            if ssm['sm'][i] != nan and ~np.isnan(ssm['sm'][i]):
                qflag[i,
                      c] = 1 + qflag[i -
                                     1][c] * np.exp(-(1. / float(t_value[c])))
                if ssm['sm_uncertainty'][i] != nan and ~np.isnan(
                        ssm['sm_uncertainty'][i]):
                    swi_noise[i, c] = sqrt(contr1_curr[c] + contr2[c] +
                                           swi_error[c])
                else:
                    swi_noise[i, c] = np.nan
            else:
                swi[i, c] = np.nan
                swi_noise[i, c] = np.nan
                qflag[i,
                      c] = qflag[i - 1][c] * np.exp(-(1. / float(t_value[c])))

    # If there are no input ssm values, gain_out should not be updated and should be the same as gain_in.
    if gain_in and len(ssm[~np.isnan(ssm['sm'])]) == 0:
        gain_out = gain_in
        gain_out['last_qflag'] = qflag[-1]  # Update to last qflag

    # If there is only 1 input SSM value:
    elif gain_in and len(ssm[~np.isnan(ssm['sm'])]) == 1:
        gain_out = {
            'last_jd': last_jd,
            'last_gain': gain_curr,
            'last_contr1': contr1_curr,
            'last_G': G_curr,
            'last_JT': JT_curr,
            'last_swi': swi[first_ssm_idx],
            'last_noise': swi_noise[first_ssm_idx],
            'last_qflag': qflag[-1]
        }  # NOTE this is the last qflag.
    else:
        gain_out = {
            'last_jd': last_jd,  # Last date where a swi calculation was done
            'last_gain': gain_curr,
            'last_contr1': contr1_curr,
            'last_G': G_curr,
            'last_JT': JT_curr,
            'last_swi': swi[last_valid_index],
            'last_noise': swi_noise[last_valid_index],
            'last_qflag': qflag[-1]
        }  # NOTE this is the last qflag.

    dtype_list = [('swi_jd', np.float64)]
    for t in t_value:
        dtype_list.append(('swi_noise_{}'.format(t), np.float32))
        dtype_list.append(('swi_{}'.format(t), np.float32))
        dtype_list.append(('qflag_{}'.format(t), np.float32))

    swi_ts = np.zeros(ssm.size, dtype=np.dtype(dtype_list))

    swi_ts['swi_jd'] = ssm['sm_jd']
    for i, t in enumerate(t_value):
        swi_ts['swi_{}'.format(t)] = swi[:, i]
        swi_ts['swi_noise_{}'.format(t)] = swi_noise[:, i]
        swi_ts['qflag_{}'.format(t)] = 100 * qflag[:, i] / qflag_norm[i]

    return swi_ts, gain_out


def calc_swi_ts(ssm_ts,
                swi_jd,
                gain_in=None,
                t_value=[1, 5, 10, 15, 20],
                nom_init=0,
                denom_init=0,
                nan=-9999.):
    """
    Time series calculation of the Soil Water Index.

    Parameters
    ----------
    ssm_ts : numpy.ndarray or dict
        Surface soil moisture time series with fields: sm_jd, sm, sm_noise
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
        gain_in = {
            'last_jd': np.float64(0),
            'denom': np.full(len(t_value), denom_init, dtype=np.float64),
            'nom': np.full(len(t_value), nom_init, dtype=np.float64),
            'nom_noise': np.zeros(len(t_value))
        }

    if gain_in['denom'][0] == 1:
        last_jd_var = ssm_ts['sm_jd'][0]
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
            swi_calc_cy_noise(ssm_ts['sm_jd'], ssm_ts['sm'], t_value, swi_jd,
                              gain_in['nom'], gain_in['denom'], last_jd_var,
                              ssm_ts['sm_noise'], nom_noise)
        swi_qflag = np.zeros_like(swi)
    else:
        swi, swi_qflag, nom, denom, last_jd_var = \
            swi_calc_cy(ssm_ts['sm_jd'], ssm_ts['sm'], t_value, swi_jd,
                        gain_in['nom'], gain_in['denom'], last_jd_var,
                        norm_factor, nan)
        swi_noise = np.zeros_like(swi)

    gain_out = {
        'denom': denom,
        'nom': nom,
        'last_jd': last_jd_var,
        'nom_noise': nom_noise
    }

    dtype_list = [('swi_jd', np.float64)]

    for t in t_value:
        dtype_list.append(('swi_{}'.format(t), np.float32))
        dtype_list.append(('swi_noise_{}'.format(t), np.float32))
        dtype_list.append(('swi_qflag_{}'.format(t), np.float32))

    swi_ts = np.zeros(swi_jd.size, dtype=np.dtype(dtype_list))

    swi_ts['swi_jd'] = swi_jd
    for i, t in enumerate(t_value):
        swi_ts['swi_{}'.format(t)] = swi[:, i]
        swi_ts['swi_noise_{}'.format(t)] = swi_noise[:, i]
        swi_ts['swi_qflag_{}'.format(t)] = swi_qflag[:, i]

    return swi_ts, gain_out


def calc_swi_noise_rec(ssm_ts, t_value, last_den=1, last_nom=0):
    """
    Recursive calculation of Soil Water Index (SWI) noise.

    Parameters
    ----------
    ssm_ts : numpy.ndarray
        Surface soil moisture time series with fields: sm_jd, sm, sm_noise
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

    last_jd = ssm_ts['sm_jd'][0]
    den_values = np.ones(len_t_value)
    nom_values = np.zeros(len_t_value)

    den_values.fill(last_den)
    nom_values.fill(last_nom)

    for i in range(0, len_sm):
        for c in range(0, len_t_value):
            tdiff = (ssm_ts['sm_jd'][i] - last_jd) / t_value[c]
            fn = exp((-1) * tdiff)
            den[i][c] = 1 + fn * den_values[c]
            den_sn = den[i][c]**2

            exp_term = exp((-2) * tdiff)
            nom_sn[i][c] = nom_values[c] * exp_term + ssm_ts['sm_noise'][i]

            swi_noise[i][c] = nom_sn[i][c] / den_sn
            den_values[c] = den[i][c]
            nom_values[c] = nom_sn[i][c]
        last_jd = ssm_ts['sm_jd'][i]

    gain_out = {
        'denom': den_values,
        'nom': nom_values,
        'last_jd': last_jd,
        'nom_noise': 0
    }

    dtype_list = [('swi_jd', np.float64)]
    for t in t_value:
        dtype_list.append(('swi_noise_{}'.format(t), np.float32))

    swi_ts = np.zeros(ssm_ts.size, dtype=np.dtype(dtype_list))

    swi_ts['swi_jd'] = ssm_ts['sm_jd']
    for i, t in enumerate(t_value):
        swi_ts['swi_noise_{}'.format(t)] = swi_noise[:, i]

    return swi_ts, gain_out
