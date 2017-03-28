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
from pyswi.swi.swi import process_swi_pd
from pyswi.swi.swi import process_swi
from pytesmo.time_series.filters import exp_filter
import pytesmo.timedate.julian as julian
import unittest
import numpy as np
import pandas as pd


class SwiTest(unittest.TestCase):

    def setUp(self):
        """
        Setup I/O objects, configuration file and grid.
        """

    def test_process_uf_pd_calc(self):
        """
        Test correct calculation of SWI, compared to the pytesmo calculation.
        """
        ctime = 5

        dates = pd.date_range('2007-01-01', periods=9).to_julian_date().values
        pd_dates = pd.date_range('2007-01-01', periods=9)

        ssm_userformat_data_test = \
            pd.DataFrame({'jd': pd_dates,
                          'sm': [10, 20, 30, 40, 50, 60, 70, 80, 90],
                          'ssf': [1, 1, 1, 1, 1, 1, 1, 1, 1]}, index=pd_dates)

        swi_ts_test, gain_test = process_swi_pd(ssm_userformat_data_test,
                                                ctime=[ctime])

        sm_in = ssm_userformat_data_test['sm'].values.astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(dates, dtype=float),
                                 ctime=ctime)

        swi_ts = swi_ts_test['SWI_005']

        np.testing.assert_array_almost_equal(swi_ts, pytesmo_swi, 4)

    def test_process_uf_pd_gain(self):
        """
        Test gain in/out of the SWI calculation.
        It calculates the first 4 days and returns the value and gain and then
        starts the calculation of day 5 to 9 and
        compares it with pytesmo over the whole period.
        """

        ctime_pytesmo = 5
        ctime = [5, 50, 100]

        dates = pd.date_range('2007-01-01', periods=9).to_julian_date().values
        pd_dates = pd.date_range('2007-01-01', periods=9)

        sm_data = [10, 20, 30, 40, 50, 60, 70, 80, 90]
        ssf_data = [1, 1, 1, 1, 1, 1, 1, 1, 1]

        ssm_ut_whole = pd.DataFrame({'jd': pd_dates, 'sm': sm_data,
                                     'ssf': ssf_data}, index=pd_dates)

        pd_dates_part1 = pd_dates[:-4]
        pd_dates_part2 = pd_dates[5:]
        sm_data_part1 = sm_data[:-4]
        sm_data_part2 = sm_data[5:]
        ssf_data_part1 = ssf_data[:-4]
        ssf_data_part2 = ssf_data[5:]

        ssm_ut_part1 = pd.DataFrame({'jd': pd_dates_part1, 'sm': sm_data_part1,
                                     'ssf': ssf_data_part1},
                                    index=pd_dates_part1)

        ssm_ut_part2 = pd.DataFrame({'jd': pd_dates_part2, 'sm': sm_data_part2,
                                     'ssf': ssf_data_part2},
                                    index=pd_dates_part2)

        swi_ts_part1, gain_test = process_swi_pd(ssm_ut_part1, ctime=ctime)

        swi_ts_part2, gain_test = process_swi_pd(ssm_ut_part2, ctime=ctime,
                                                 gain_in=gain_test)

        sm_in = ssm_ut_whole['sm'].values.astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(dates, dtype=float),
                                 ctime=ctime_pytesmo)

        swi_ts_ctime1 = swi_ts_part1['SWI_005'].tolist()
        swi_ts_ctime2 = swi_ts_part2['SWI_005'].tolist()

        np.testing.assert_array_almost_equal(swi_ts_ctime1+swi_ts_ctime2,
                                             pytesmo_swi, 5)

    def test_process_uf_pd_time(self):
        """
        Test correct calculation of SWI, compared to the pytesmo calculation.
        """
        ctime = 5

        pd_dates = pd.date_range('2007-01-01', periods=9)

        dates_area = \
            pd.date_range('2007-01-03', periods=3).to_julian_date().values

        date_from = min(dates_area)
        date_to = max(dates_area)

        proc_param = {'date_from': date_from, 'date_to': date_to}

        ssm_userformat_data_test = \
            pd.DataFrame({'jd': pd_dates,
                          'sm': [10, 20, 30, 40, 50, 60, 70, 80, 90],
                          'ssf': [1, 1, 1, 1, 1, 1, 1, 1, 1]}, index=pd_dates)

        sm_pytesmo = np.array([30, 40, 50])

        swi_ts_test, gain_test = \
            process_swi_pd(ssm_userformat_data_test, ctime=[ctime],
                           proc_param=proc_param)

        sm_in = sm_pytesmo.astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(dates_area, dtype=float),
                                 ctime=ctime)

        swi_ts = swi_ts_test['SWI_005']

        np.testing.assert_array_almost_equal(swi_ts, pytesmo_swi, 4)

    def test_process_swi_calc(self):
        """
        Test correct calculation of SWI, compared to the pytesmo calculation.
        """
        ctime_pytesmo = 5
        ctime = [5, 50, 100]

        dates = pd.date_range('2007-01-01', periods=9).to_julian_date().values

        sm = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])

        swi_ts_test, gain_test = process_swi(sm, dates, ctime=ctime)

        sm_in = sm.astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(dates, dtype=float),
                                 ctime=ctime_pytesmo)

        swi_ts = swi_ts_test['SWI_005']

        np.testing.assert_array_almost_equal(swi_ts, pytesmo_swi, 4)

    def test_process_swi_gain(self):
        """
        Test gain in/out of the SWI calculation.
        It calculates the first 4 days and returns the value and gain and
        then starts the calculation of day 5 to 9 and
        compares it with pytesmo over the whole period.
        """
        ctime = 5

        dates = pd.date_range('2007-01-01', periods=9).to_julian_date().values

        sm_data = [10, 20, 30, 40, 50, 60, 70, 80, 90]

        dates_part1 = dates[:-4]
        dates_part2 = dates[5:]
        sm_data_part1 = np.array(sm_data[:-4])
        sm_data_part2 = np.array(sm_data[5:])

        swi_ts_part1, gain_test = process_swi(sm_data_part1, dates_part1,
                                              ctime=[ctime])

        swi_ts_part2, gain_test = process_swi(sm_data_part2, dates_part2,
                                              ctime=[ctime], gain_in=gain_test)

        sm_in = np.array(sm_data).astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(dates, dtype=float),
                                 ctime=ctime)

        swi_ts_ctime1 = swi_ts_part1['SWI_005'].tolist()
        swi_ts_ctime2 = swi_ts_part2['SWI_005'].tolist()

        np.testing.assert_array_almost_equal(swi_ts_ctime1+swi_ts_ctime2,
                                             pytesmo_swi, 5)

    def test_process_swi_time(self):
        """
        Test correct calculation of SWI, compared to the pytesmo calculation.
        """
        ctime = 5

        dates = pd.date_range('2007-01-01', periods=9).to_julian_date().values

        dates_area = \
            pd.date_range('2007-01-03', periods=3).to_julian_date().values

        date_from = min(dates_area)
        date_to = max(dates_area)

        proc_param = {'date_from': date_from, 'date_to': date_to}

        sm = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])

        sm_pytesmo = np.array([30, 40, 50])

        swi_ts_test, gain_test = process_swi(sm, dates, ctime=[ctime],
                                             proc_param=proc_param)

        sm_in = sm_pytesmo.astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(dates_area, dtype=float),
                                 ctime=ctime)

        swi_ts = swi_ts_test['SWI_005']

        np.testing.assert_array_almost_equal(swi_ts, pytesmo_swi, 4)

    def test_process_swi_not_daily_out(self):
        """
        Test correct calculation of SWI, compared to the pytesmo calculation.
        """
        ctime_pytesmo = 5
        ctime = [5, 50, 100]

        dates = np.array([2454102.3093969999, 2454103.3644969999,
                          2454103.873199, 2454104.9282769999,
                          2454105.9139760002, 2454106.3213539999,
                          2454107.3070530002, 2454108.3621530002,
                          2454108.8708330002, 2454109.9259549999])

        sm = np.array([84, 88, 79, 87, 91, 93, 92, 85, 85, 90])

        swi_ts_test, gain_test = process_swi(sm, dates, ctime=ctime)

        sm_in = sm.astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(dates, dtype=float),
                                 ctime=ctime_pytesmo)

        swi_ts = swi_ts_test['SWI_005']

        np.testing.assert_array_almost_equal(swi_ts, pytesmo_swi, 4)

    def test_process_swi_daily_out(self):
        """
        Test correct calculation of SWI, compared to the pytesmo calculation.
        """
        ctime_pytesmo = 5
        ctime = [5]

        dates = np.array([2454102.3093969999, 2454103.3644969999,
                          2454103.873199, 2454104.9282769999,
                          2454105.9139760002, 2454106.3213539999,
                          2454107.3070530002, 2454108.3621530002,
                          2454108.8708330002, 2454109.9259549999])

        pytesmo_dates = np.array([2454102.3093969999, 2454103.3644969999,
                                  2454103.873199, 2454104.9282769999,
                                  2454105.9139760002, 2454106.3213539999,
                                  2454107.3070530002, 2454108.3621530002,
                                  2454108.8708330002, 2454109.9259549999])

        sm = np.array([84, 88, 79, 87, 91, 93, 92, 85, 85, 90])

        swi_ts_test, gain_test = process_swi(sm, dates, ctime=ctime,
                                             jd_daily_out=True)

        pytesmo_sm = np.array([84, 88, 79, 87, 91, 93, 92, 85, 85, 90])

        sm_in = pytesmo_sm.astype(float)

        pytesmo_swi = exp_filter(sm_in, np.asarray(pytesmo_dates, dtype=float),
                                 ctime=ctime_pytesmo)

        swi_ts = swi_ts_test['SWI_005']

        months, days, years = julian.caldat(swi_ts_test['jd'])

        np.testing.assert_array_equal(years, [2007, 2007, 2007, 2007, 2007,
                                              2007, 2007, 2007])
        np.testing.assert_equal(days, range(1, 9))
        np.testing.assert_equal(months, [1, 1, 1, 1, 1, 1, 1, 1])

        pytesmo_swi = np.delete(pytesmo_swi, [4, 9])

        np.testing.assert_array_almost_equal(swi_ts, pytesmo_swi, 4)

if __name__ == '__main__':
    unittest.main()
