# Copyright (c) 2020, Vienna University of Technology (TU Wien), Department
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

from pyswi.swi_img.iterative_swi import IterativeMultiSWI, IterativeSWI, calc_swi


class IterativeSwiTest(unittest.TestCase):

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

        next_ssm = np.array([20., 20., 20.])
        next_ssm_jd = np.array([1, 2, 3])
        tvalue = 5
        swi = np.array([10, 10, 10])
        gain = np.array([1., 1., 1.])
        jd = np.array([0, 0, 0])

        swi_ref = np.array([15.49833997, 15.9868766, 16.45656306])
        gain_ref = np.array([0.549834, 0.59868766, 0.64565631])

        swi, gain = calc_swi(next_ssm, next_ssm_jd, tvalue, swi,
                                         jd, gain)

        np.testing.assert_array_almost_equal(swi, swi_ref, decimal=6)
        np.testing.assert_array_almost_equal(gain, gain_ref, decimal=6)

    def test_iterative_swi(self):
        """
        Test correct calculation of SWI including the import of the test data
        compared to a hardcoded calculated output.
        """

        ssm = np.array([45., 46., 32., 51., 65., 23., 54., 44., 23., 42., 46.,
                        75., 43., 23., 12., 23., 11., 9., 23., 24.])
        ssm_jd = np.array(
            [2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456598.64583333, 2456598.64583333, 2456599.64583333,
             2456599.64583333, 2456599.64583333])
        tvalue = 5

        swi_ref = np.array([24.7425298791, 25.2923638764, 17.594687914,
                            28.0415338629, 35.7392098253, 12.6461819382,
                            29.6910358549, 24.1926958817, 12.6461819382,
                            23.0930278871, 25.2923638764, 41.2375497984,
                            23.6428618844, 12.6461819382, 6.5980079677,
                            13.7698161826, 6.5855642612, 5.810906756,
                            14.8500950432, 15.4957513494])

        iter_swi = IterativeSWI(ssm, os.path.join(self.curpath(),
                                                  'data',
                                                  'iterative_swi'), tvalue)

        swi = iter_swi.calc_iter(ssm_jd, ssm)

        np.testing.assert_array_almost_equal(swi, swi_ref, decimal=6)

    def test_iterative_swi_load_save(self):
        """
        Test correct save and load process of the iterative SWI.
        """
        out_path = mkdtemp()

        ssm = np.array([45., 46., 32., 51., 65., 23., 54., 44., 23., 42., 46.,
                        75., 43., 23., 12., 23., 11., 9., 23., 24.])

        ssm_jd = np.array(
            [2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456598.64583333, 2456598.64583333, 2456599.64583333,
             2456599.64583333, 2456599.64583333])

        tvalue = 5

        swi_ref = np.array([24.7425298791, 25.2923638764, 17.594687914,
                            28.0415338629, 35.7392098253, 12.6461819382,
                            29.6910358549, 24.1926958817, 12.6461819382,
                            23.0930278871, 25.2923638764, 41.2375497984,
                            23.6428618844, 12.6461819382, 6.5980079677,
                            13.7698161826, 6.5855642612, 5.810906756,
                            14.8500950432, 15.4957513494])

        iter_swi = IterativeSWI(ssm, os.path.join(self.curpath(),
                                                  'data',
                                                  'iterative_swi'), tvalue)

        swi = iter_swi.calc_iter(ssm_jd, ssm)

        iter_swi.iterstepdata.path = out_path

        iter_swi.store_iter_data()

        np.testing.assert_array_almost_equal(swi, swi_ref, decimal=6)

        iter_swi2 = IterativeSWI(ssm, out_path, tvalue)

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

        swi2 = iter_swi2.calc_iter(ssm_jd2, ssm2)

        np.testing.assert_array_almost_equal(swi2, swi_ref2, decimal=6)

    def test_multi_iterative_swi(self):
        """
        Test correct calculation of SWI including the import of the test data
        compared to a hardcoded calculated output for different tvalue values.
        """

        ssm = np.array([45., 46., 32., 51., 65., 23., 54., 44., 23., 42., 46.,
                        75., 43., 23., 12., 23., 11., 9., 23., 24.])
        ssm_jd = np.array(
            [2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456598.64583333, 2456598.64583333, 2456599.64583333,
             2456599.64583333, 2456599.64583333])
        tvalues = np.array([5, 50, 75])

        swi_005_ref = np.array([24.7425298791, 25.2923638764, 17.594687914,
                                28.0415338629, 35.7392098253, 12.6461819382,
                                29.6910358549, 24.1926958817, 12.6461819382,
                                23.0930278871, 25.2923638764, 41.2375497984,
                                23.6428618844, 12.6461819382, 6.5980079677,
                                13.7698161826, 6.5855642612, 5.810906756,
                                14.8500950432, 15.4957513494])

        swi_050_ref = np.array([22.7249925003, 23.2299923336, 16.1599946669,
                                25.7549915003, 32.8249891671, 11.6149961668,
                                27.2699910004, 22.219992667, 11.6149961668,
                                21.2099930003, 23.2299923336, 37.8749875005,
                                21.7149928336, 11.6149961668, 6.0599980001,
                                11.7299693382, 5.6099853357, 4.6349595146,
                                11.8448965372, 12.3598920389])

        swi_075_ref = np.array([22.6499977778, 23.1533310618, 16.1066650864,
                                25.6699974815, 32.7166634568, 11.5766655309,
                                27.1799973334, 22.1466644939, 11.5766655309,
                                21.139997926, 23.1533310618, 37.7499962964,
                                21.6433312099, 11.5766655309, 6.0399994074,
                                11.6533242476, 5.573328988, 4.5899880019,
                                11.7299693382, 12.2399680051])

        iter_swi = IterativeMultiSWI(ssm, os.path.join(self.curpath(),
                                     'data', 'iterative_swi'), tvalues)

        swi = iter_swi.calc_iter(ssm_jd, ssm)

        np.testing.assert_array_almost_equal(swi['SWI_005'],
                                             swi_005_ref, decimal=6)
        np.testing.assert_array_almost_equal(swi['SWI_050'],
                                             swi_050_ref, decimal=6)
        np.testing.assert_array_almost_equal(swi['SWI_075'],
                                             swi_075_ref, decimal=6)

    def test_iterative_multi_swi_load_save(self):
        """
        Test correct save and load process of the iterative SWI.
        """
        out_path = mkdtemp()

        ssm = np.array([45., 46., 32., 51., 65., 23., 54., 44., 23., 42., 46.,
                        75., 43., 23., 12., 23., 11., 9., 23., 24.])

        ssm_jd = np.array(
            [2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456597.64583333, 2456597.64583333, 2456597.64583333,
             2456598.64583333, 2456598.64583333, 2456599.64583333,
             2456599.64583333, 2456599.64583333])
        tvalues = np.array([5, 50, 75])

        iter_swi = IterativeMultiSWI(ssm, os.path.join(self.curpath(),
                                     'data', 'iterative_swi'), tvalues)

        iter_swi.calc_iter(ssm_jd, ssm)

        for iter in iter_swi.calc_swi:
            iter.iterstepdata.path = out_path

        iter_swi.store_iter_data()

        iter_swi2 = IterativeMultiSWI(ssm, out_path, tvalues)

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

        swi_005_ref = np.array([32.881162537, 33.6118550378, 23.3821600263,
                                37.2653175419, 47.4950125534, 16.8059275189,
                                39.4573950444, 32.1504700362, 16.8059275189,
                                30.6890850345, 33.6118550378, 54.8019375617,
                                31.4197775354, 16.8059275189, 8.7683100099,
                                17.668451265, 8.4501288658, 7.2169955022,
                                18.44343295, 19.2453213392])

        swi_050_ref = np.array([30.2989801017, 30.9722907707, 21.5459414057,
                                34.3388441153, 43.7651934803, 15.4861453853,
                                36.3587761221, 29.6256694328, 15.4861453853,
                                28.2790480949, 30.9722907707, 50.4983001695,
                                28.9523587639, 15.4861453853, 8.0797280271,
                                15.5869756945, 7.4546405495, 6.1384272595,
                                15.6870918854, 16.3691393586])

        swi_075_ref = np.array([30.1995496496, 30.8706507529, 21.4752353064,
                                34.2261562696, 43.6215717161, 15.4353253765,
                                36.2394595795, 29.5284485463, 15.4353253765,
                                28.1862463396, 30.8706507529, 50.3325827494,
                                28.857347443, 15.4353253765, 8.0532132399,
                                15.5028591278, 7.4144108872, 6.0926400726,
                                15.5700801855, 16.2470401935])

        swi2 = iter_swi2.calc_iter(ssm_jd2, ssm2)

        np.testing.assert_array_almost_equal(swi2['SWI_005'],
                                             swi_005_ref, decimal=6)
        np.testing.assert_array_almost_equal(swi2['SWI_050'],
                                             swi_050_ref, decimal=6)
        np.testing.assert_array_almost_equal(swi2['SWI_075'],
                                             swi_075_ref, decimal=6)

if __name__ == '__main__':
    unittest.main()