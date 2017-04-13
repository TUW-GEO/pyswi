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
from pyswi.swi.swi import iterative_swi
from pyswi.swi.iterative_swi import IterativeSWI
import unittest
import numpy as np
import os
import pandas as pd
from tempfile import mkdtemp

class SwiTest(unittest.TestCase):

    def curpath(self):
        return os.path.abspath(os.path.dirname(__file__))

    def setUp(self):
        """
        Setup I/O objects, configuration file and grid.
        """

    def test_iterative_swi_calc(self):
        """
        Test correct calculation of SWI to a hardcoded calculated output.
        """

        ssm = np.array([20., 20., 20.])
        ssm_jd = np.array([1, 2, 3])
        time_scale = 5
        prev_swi = np.array([10, 10, 10])
        prev_gain = np.array([1., 1., 1.])
        prev_qflag = np.array([40., 40., 40.])
        prev_jd = np.array([0, 0, 0])

        swi_ref = np.array([15.49833997, 15.9868766, 16.45656306])
        qflag_ref = np.array([33.74923012, 27.81280184, 22.95246544])
        gain_ref = np.array([0.549834, 0.59868766, 0.64565631])

        swi, qflag, gain = iterative_swi(ssm, ssm_jd, time_scale, prev_swi,
                                             prev_gain, prev_qflag, prev_jd)

        np.testing.assert_array_almost_equal(swi, swi_ref, decimal=6)
        np.testing.assert_array_almost_equal(qflag, qflag_ref, decimal=6)
        np.testing.assert_array_almost_equal(gain, gain_ref, decimal=6)

    def test_iterative_swi_calc(self):
        """
        Test correct calculation of SWI to a hardcoded calculated output.
        """

        ssm = np.array([45., 46., 32., 51., 65., 23., 54., 44., 23., 42., 46.,
                        75., 43., 23., 12., 23., 11., 9., 23., 24.])
        ssm_jd = np.array([2456597.64583333, 2456597.64583333, 2456597.64583333,
                           2456597.64583333, 2456597.64583333, 2456597.64583333,
                           2456597.64583333, 2456597.64583333, 2456597.64583333,
                           2456597.64583333, 2456597.64583333, 2456597.64583333,
                           2456597.64583333, 2456597.64583333, 2456597.64583333,
                           2456598.64583333, 2456598.64583333, 2456599.64583333,
                           2456599.64583333, 2456599.64583333])
        ctime = 5

        swi_ref = np.array([24.7425298791, 25.2923638764, 17.594687914,
                            28.0415338629, 35.7392098253, 12.6461819382,
                            29.6910358549, 24.1926958817, 12.6461819382,
                            23.0930278871, 25.2923638764, 41.2375497984,
                            23.6428618844, 12.6461819382, 6.5980079677,
                            13.7698161826, 6.5855642612, 5.810906756,
                            14.8500950432, 15.4957513494])

        # gain_ref = np.array([0.5498339973, 0.5498339973, 0.5498339973,
        #                      0.5498339973, 0.5498339973, 0.5498339973,
        #                      0.5498339973, 0.5498339973, 0.5498339973,
        #                      0.5498339973, 0.5498339973, 0.5498339973,
        #                      0.5498339973, 0.5498339973, 0.5498339973,
        #                      0.5986876601, 0.5986876601, 0.6456563062,
        #                      0.6456563062, 0.6456563062])

        outdir = os.path.join(mkdtemp())
        print outdir
        iter_swi = IterativeSWI(ssm, os.path.join(self.curpath(),
                                         'data',
                                         'iterative_swi'))

        swi = iter_swi.calc_iter(ssm_jd, ssm, ctime)

        np.testing.assert_array_almost_equal(swi, swi_ref, decimal=6)

    # def test_iterative_swi_process(self):
    #     """
    #     Test correct calculation of SWI, compared to the pytesmo calculation.
    #     """
    #     outdir = os.path.join(mkdtemp())
    #     print outdir
    #     ssm = np.array([20., 20., 20.])
    #     ssm_jd = pd.date_range('2007-01-01', periods=3).to_julian_date().values
    #
    #     ctime = 5
    #
    #     iter_swi = IterativeSWI(ssm, outdir)
    #
    #     swi = iter_swi.calc_iter(ssm_jd, ssm, ctime)
    #
    #     swi_ref = np.array([15.49833997, 15.9868766, 16.45656306])
    #
    #     iter_swi.store_iter_data()
    #
    #     iter_swi.load_iter_data()
    #
    #     ssm = np.array([10., 10., 25.])
    #     ssm_jd = pd.date_range('2007-01-03', periods=3).to_julian_date().values
    #
    #     swi = iter_swi.calc_iter(ssm_jd, ssm, ctime)
    #
    #     np.testing.assert_array_almost_equal(swi, swi_ref, decimal=6)


if __name__ == '__main__':
    unittest.main()
