# Copyright (c) 2019, Vienna University of Technology (TU Wien), Department
# of Geodesy and Geoinformation (GEO).
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

"""
This module tests the soil water index calculation from SSM time series.
"""

import numpy as np
from numpy.lib.recfunctions import unstructured_to_structured
import pandas as pd

from pyswi.swi_ts.swi_ts import calc_swi_ts, calc_swi_noise_rec


def test_process_swi_calc():
    """
    Test correct calculation of SWI, compared to a hardcoded calculation output.
    """
    t_value = [5, 50, 100]

    swi_jd = pd.date_range(
        '2007-01-01', periods=9).to_julian_date().values.astype(np.float64)

    sm = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90], dtype=np.float32)

    dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32)])
    ssm_ts = unstructured_to_structured(
        np.hstack((swi_jd[:, np.newaxis], sm[:, np.newaxis])), dtype=dtype)

    swi_ts, gain_out = calc_swi_ts(ssm_ts, swi_jd, t_value=t_value)

    swi_ref = {
        'swi_5': np.array([10., 15.49834, 21.32452, 27.472094, 33.93228, 40.69421,
                           47.7452, 55.07107, 62.65647]),
        'swi_50': np.array([10., 15.049998, 20.133324, 25.249971, 30.399931, 35.58319,
                            40.799732, 46.049545, 51.332603]),
        'swi_100': np.array([10., 15.025, 20.066666, 25.124996, 30.199991, 35.29165,
                             40.399967, 45.524944, 50.666576])
    }

    for t in t_value:
        np.testing.assert_array_almost_equal(swi_ts['swi_{}'.format(t)],
                                             swi_ref['swi_{}'.format(t)], 4)

def test_process_swi_calc_nan_values():
    """
    Test correct calculation of SWI, compared to a hardcoded calculation output
    including nan values.
    """
    t_value = [5, 50, 100]

    swi_jd = pd.date_range(
        '2007-01-01', periods=9).to_julian_date().values.astype(np.float64)

    sm = np.array([10, 20, 30, 40, 50, 60, 999, 80, 10])

    dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32)])
    ssm_ts = unstructured_to_structured(
        np.hstack((swi_jd[:, np.newaxis], sm[:, np.newaxis])), dtype=dtype)

    swi_ts, gain_out = calc_swi_ts(ssm_ts, swi_jd, t_value=t_value, nan=999)

    swi_ref = {
        'swi_5': np.array([10., 15.49834, 21.32452, 27.472094, 33.93228, 40.69421,
                           np.nan, 51.660824, 41.072067]),
        'swi_50': np.array([10., 15.049998, 20.133324, 25.249971, 30.399931, 35.58319,
                            np.nan, 42.430466, 38.023155]),
        'swi_100': ([10., 15.025, 20.066666, 25.124996, 30.199991, 35.29165,
                     np.nan, 41.928066, 37.76523])
    }

    for t in t_value:
        np.testing.assert_array_almost_equal(swi_ts['swi_{}'.format(t)],
                                             swi_ref['swi_{}'.format(t)], 4)

def test_process_swi_gain():
    """
    Test gain in/out of the SWI calculation.
    It calculates the first 4 days and returns the value and gain and
    then starts the calculation of day 5 to 9 and
    compares it with hardcoded output over the whole period.
    """
    t_value = [5, 50, 100]

    swi_jd = pd.date_range(
        '2007-01-01', periods=9).to_julian_date().values.astype(np.float64)

    sm = np.array([10, 20, 30, 40, 50, 60, 999, 80, 10])

    dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32)])
    ssm_ts_part1 = unstructured_to_structured(
        np.hstack((swi_jd[:5, np.newaxis], sm[:5, np.newaxis])), dtype=dtype)
    ssm_ts_part2 = unstructured_to_structured(
        np.hstack((swi_jd[5:, np.newaxis], sm[5:, np.newaxis])), dtype=dtype)

    swi_ts_1, gain_1 = calc_swi_ts(ssm_ts_part1, swi_jd[:5],
                                   t_value=t_value, nan=999)

    swi_ts_2, gain_2 = calc_swi_ts(ssm_ts_part2, swi_jd[5:],
                                   t_value=t_value, gain_in=gain_1,
                                   nan=999)

    swi_ref = {
        'swi_5': np.array([10., 15.49834, 21.32452, 27.472094, 33.93228, 40.69421,
                           np.nan, 51.660824, 41.072067]),
        'swi_50': np.array([10., 15.049998, 20.133324, 25.249971, 30.399931, 35.58319,
                            np.nan, 42.430466, 38.023155]),
        'swi_100': ([10., 15.025, 20.066666, 25.124996, 30.199991, 35.29165,
                     np.nan, 41.928066, 37.76523])
    }

    for t in t_value:
        comb_swi_ts = swi_ts_1['swi_{}'.format(t)].tolist() + \
                      swi_ts_2['swi_{}'.format(t)].tolist()
        np.testing.assert_array_almost_equal(comb_swi_ts, swi_ref['swi_{}'.format(t)], 4)

def test_process_swi_not_daily_out():
    """
    Test correct calculation of SWI with irregular timestamps,
    compared to the a hardcoded calculation output.
    """
    t_value = [5, 50, 100]

    swi_jd = np.array([2454102.3093969999, 2454103.3644969999,
                       2454103.873199, 2454104.9282769999,
                       2454105.9139760002, 2454106.3213539999,
                       2454107.3070530002, 2454108.3621530002,
                       2454108.8708330002, 2454109.9259549999])

    sm = np.array([84, 88, 79, 87, 91, 93, 92, 85, 85, 90])

    dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32)])
    ssm_ts = unstructured_to_structured(
        np.hstack((swi_jd[:, np.newaxis], sm[:, np.newaxis])), dtype=dtype)

    swi_ts, gain_out = calc_swi_ts(ssm_ts, swi_jd, t_value=t_value)

    swi_ref = {
        'swi_5': np.array([84., 86.21024, 83.47359, 84.59898, 86.39057, 87.93006,
                           88.829475, 88.008446, 87.43135, 87.9233]),
        'swi_50': np.array([84., 86.0211, 83.64838, 84.50836, 85.854836, 87.090576,
                            87.826706, 87.45131, 87.16135, 87.46732]),
        'swi_100': np.array([84., 86.01055, 83.65755, 84.50414, 85.827286, 87.04517,
                             87.77045, 87.41343, 87.13668, 87.434044])
    }

    for t in t_value:
        np.testing.assert_array_almost_equal(swi_ts['swi_{}'.format(t)],
                                             swi_ref['swi_{}'.format(t)], 4)

def test_process_swi_noise():
    """
    Test correct calculation of SWI Noise inside the SWI calculation
    process.
    """
    t_value = [5, 50, 100]

    swi_jd = pd.date_range(
        '2007-01-01', periods=9).to_julian_date().values.astype(np.float64)

    sm_noise = np.zeros(9)
    sm_noise.fill(3)

    sm = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])

    dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32),
                      ('sm_noise', np.float32)])
    ssm_ts = unstructured_to_structured(
        np.hstack((swi_jd[:, np.newaxis], sm[:, np.newaxis],
                   sm_noise[:, np.newaxis])), dtype=dtype)

    swi_ts1, _ = calc_swi_ts(ssm_ts, swi_jd, t_value=t_value,
                             nom_init=1, denom_init=1)

    swi_ts2, _ = calc_swi_noise_rec(ssm_ts, t_value)

    for t in t_value:
        np.testing.assert_array_almost_equal(
            swi_ts1['swi_noise_{}'.format(t)],
            swi_ts2['swi_noise_{}'.format(t)], 4)

def test_swi_gain_noise():
    """
    Test correct calculation of SWI noise gain.
    """
    t_value = [5, 50, 100]

    swi_jd = pd.date_range('2007-01-01', periods=9).to_julian_date().values
    sm = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])
    sm_noise = np.zeros(9)
    sm_noise.fill(3)

    dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32),
                      ('sm_noise', np.float32)])
    ssm_ts = unstructured_to_structured(
        np.hstack((swi_jd[:, np.newaxis], sm[:, np.newaxis],
                   sm_noise[:, np.newaxis])), dtype=dtype)

    swi_ts_all, gain_all = calc_swi_ts(ssm_ts, swi_jd, t_value=t_value)
    swi_ts_1, gain_1 = calc_swi_ts(ssm_ts[:5], swi_jd[:5], t_value=t_value)
    swi_ts_2, gain_2 = calc_swi_ts(ssm_ts[5:], swi_jd[5:],
                                   t_value=t_value, gain_in=gain_1)

    for t in t_value:
        comb_swi_ts = swi_ts_1['swi_noise_{}'.format(t)].tolist() + \
            swi_ts_2['swi_noise_{}'.format(t)].tolist()
        np.testing.assert_array_almost_equal(
            comb_swi_ts, swi_ts_all['swi_noise_{}'.format(t)], 4)
