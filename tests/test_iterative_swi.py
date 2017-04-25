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
import os
import unittest
from tempfile import mkdtemp

import numpy as np

from pyswi.swi_img.iterative_swi import IterativeSWI
from pyswi.swi_img.swi_img import iterative_swi


class SwiTest(unittest.TestCase):

    def curpath(self):
        return os.path.abspath(os.path.dirname(__file__))

    def setUp(self):
        """
        Setup I/O objects, configuration file and grid.
        """

    def test_iterative_swi_calc(self):
        """
        Test correct calculation of SWI compared to a hardcoded calculated
        output.
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

    def test_iterative_swi(self):
        """
        Test correct calculation of SWI including the import of the test data
        compared to a hardcoded calculated output.
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

        iter_swi = IterativeSWI(ssm, os.path.join(self.curpath(),
                                                  'data',
                                                  'iterative_swi'))

        swi = iter_swi.calc_iter(ssm_jd, ssm, ctime)

        np.testing.assert_array_almost_equal(swi, swi_ref, decimal=6)

    def test_iterative_swi_load_save(self):
        """
        Test correct save and load process of the iterative SWI.
        """
        out_path = mkdtemp()
        print out_path

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

        iter_swi = IterativeSWI(ssm, os.path.join(self.curpath(),
                                                  'data',
                                                  'iterative_swi'))

        swi = iter_swi.calc_iter(ssm_jd, ssm, ctime)

        iter_swi.iterstepdata.path = out_path

        iter_swi.store_iter_data()

        np.testing.assert_array_almost_equal(swi, swi_ref, decimal=6)

        iter_swi2 = IterativeSWI(ssm, out_path)

        ssm2 = np.array([45., 46., 32., 51., 65., 23., 54., 44., 23., 42., 46.,
                        75., 43., 23., 12., 23., 11., 9., 23., 24.])

        ssm_jd2 = np.array([2456598.64583333, 2456598.64583333,
                            2456598.64583333, 2456598.64583333,
                            2456598.64583333, 2456598.64583333,
                            2456598.64583333, 2456598.64583333,
                            2456598.64583333, 2456598.64583333,
                            2456598.64583333, 2456598.64583333,
                            2456598.64583333, 2456598.64583333,
                            2456598.64583333, 2456599.64583333,
                            2456599.64583333, 2456600.64583333,
                            2456600.64583333, 2456600.64583333])

        swi_ref2 = np.array([32.881162537, 33.6118550378, 23.3821600263,
                             37.2653175419, 47.4950125534, 16.8059275189,
                             39.4573950444, 32.1504700362, 16.8059275189,
                             30.6890850345, 33.6118550378, 54.8019375617,
                             31.4197775354, 16.8059275189, 8.7683100099,
                             17.668451265, 8.4501288658, 7.2169955022,
                             18.44343295, 19.2453213392])

        swi2 = iter_swi2.calc_iter(ssm_jd2, ssm2, ctime)

        np.testing.assert_array_almost_equal(swi2, swi_ref2, decimal=6)


if __name__ == '__main__':
    unittest.main()
