# Copyright (c) 2017, Vienna University of Technology (TU Wien), Department
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
This module tests the soil water index calculation.
"""
from pyswi.swi.swi import process
from pytesmo.time_series.filters import exp_filter
import unittest
import numpy as np
import pandas as pd
import datetime

class SwiTest(unittest.TestCase):

    def setUp(self):
        """
        Setup I/O objects, configuration file and grid.
        """

    def test_swi_calc(self):
        """
        Test correct calculation of SWI, compared to the pytesmo calculation.
        """
        ctime = 5

        proc_param = {}

        swi_param = {'ctime': str([ctime])}

        dates = pd.date_range('2007-01-01', periods=9).to_julian_date().values
        pd_dates = pd.date_range('2007-01-01', periods=9)

        ssm_userformat_data_test = pd.DataFrame({'jd': pd_dates, 'sm': [10, 20, 30, 40, 50, 60, 70, 80, 90],
                                                 'ssf': [1, 1, 1, 1, 1, 1, 1, 1, 1]}, index=pd_dates)

        swi_ts_test, gain_test = process(ssm_userformat_data_test, proc_param, swi_param)

        sm_in = ssm_userformat_data_test['sm'].values.astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(dates, dtype=float), ctime=ctime)

        swi_ts = swi_ts_test['SWI_005']

        np.testing.assert_array_almost_equal(swi_ts, pytesmo_swi, 4)


    def test_gain(self):
        """
        Test gain in/out of the SWI calculation.
        It calculates the first 4 days and returns the value and gain and then starts the calculation of day 5 to 9 and
        compares it with pytesmo over the whole period.
        """
        ctime = 5

        proc_param = {}

        swi_param = {'ctime': str([ctime])}

        dates = pd.date_range('2007-01-01', periods=9).to_julian_date().values
        pd_dates = pd.date_range('2007-01-01', periods=9)

        sm_data = [10, 20, 30, 40, 50, 60, 70, 80, 90]
        ssf_data = [1, 1, 1, 1, 1, 1, 1, 1, 1]

        ssm_ut_whole = pd.DataFrame({'jd': pd_dates, 'sm': sm_data, 'ssf': ssf_data}, index=pd_dates)

        pd_dates_part1 = pd_dates[:-4]
        pd_dates_part2 = pd_dates[5:]
        sm_data_part1 = sm_data[:-4]
        sm_data_part2 = sm_data[5:]
        ssf_data_part1 = ssf_data[:-4]
        ssf_data_part2 = ssf_data[5:]

        ssm_ut_part1 = pd.DataFrame({'jd': pd_dates_part1, 'sm': sm_data_part1,
                                     'ssf': ssf_data_part1}, index=pd_dates_part1)

        ssm_ut_part2 = pd.DataFrame({'jd': pd_dates_part2, 'sm': sm_data_part2,
                                     'ssf': ssf_data_part2}, index=pd_dates_part2)

        swi_ts_part1, gain_test = process(ssm_ut_part1, proc_param, swi_param)

        swi_ts_part2, gain_test = process(ssm_ut_part2, proc_param, swi_param, gain_test)

        sm_in = ssm_ut_whole['sm'].values.astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(dates, dtype=float), ctime=ctime)

        swi_ts_ctime1 = swi_ts_part1['SWI_005'].tolist()
        swi_ts_ctime2 = swi_ts_part2['SWI_005'].tolist()

        np.testing.assert_array_almost_equal(swi_ts_ctime1+swi_ts_ctime2, pytesmo_swi, 5)

    def test_ssf_frozen(self):
        """
        Test correct masking of frozen data for the SWI calculation.
        """
        ctime = 5

        proc_param = {}

        swi_param = {'ctime': str([ctime])}

        dates = pd.date_range('2007-01-01', periods=6).to_julian_date().values
        pd_dates = pd.date_range('2007-01-03', periods=9)

        ssm_userformat_data_test = pd.DataFrame({'jd': pd_dates, 'sm': [10, 20, 30, 40, 50, 60, 70, 80, 90],
                                                 'ssf': [0, 2, 3, 1, 1, 1, 1, 1, 1]}, index=pd_dates)

        swi_ts_test, gain_test = process(ssm_userformat_data_test, proc_param, swi_param)

        sm_in = np.array([40, 50, 60, 70, 80, 90])
        sm_in = sm_in.astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(dates, dtype=float), ctime=ctime)

        swi_ts_ctime5 = swi_ts_test['SWI_005']

        np.testing.assert_array_almost_equal(swi_ts_ctime5[3:], pytesmo_swi, 4)


if __name__ == '__main__':
    unittest.main()
