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
import os
from datetime import datetime

import numpy as np
from pyswi.swi_img.iterative_storage import IterStepData
from cadati.jd_date import julian2datetime

def calc_swi(next_ssm, next_ssm_jd, tvalue, swi, jd,
                  gain=None, weighted=False, den=None, n=None, wsum=None, next_w=None):

    """
    Takes input data and previous data and calculates swi_img
    values for the next date. All arrays have to be the same shape
    so that the whole batch can be calculated at once
    initializing and masking has to be done outside of this basic
    function

    Parameters
    ----------
    next_ssm : numpy.array
        surface soil moisture
    next_ssm_jd : numpy.array
        observation time of ssm given as julian date
    tvalue : int
        T-value for SWI calculation
    swi : numpy.array
        values of SWI at previous iteration
    jd : numpy.array
        values of julian date from previous iteration
    gain : numpy.array (optional, default=None)
        values of gain from previous iteration. Required for non-weighted swi calculation.
    weighted : boolean (default=False)
        Weighted swi calculation following equation (9) in the paper:
        Bauer-Marschallinger B, Paulik C, Hochstoeger S, Mistelbauer T,
        Modanesi S, Ciabatta L, Massari C, Brocca L, Wagner W.
        "Soil moisture from fusion of scatterometer and SAR:
        Closing the scale gap with temporal filtering."
        Remote Sensing. 2018 Jul;10(7):1030.
    den: numpy.array (optional, default=None)
        values of denominator from previous iteration. Required for weighted swi calculation.
    n : numpy.array (optional, default=None)
        values of total number of previous iterations. Required for weighted swi calculation.
    wsum : numpy.array (optional, default=None)
        values of sum of weights applied in previous iteration. Required for weighted swi calculation.
    next_w : numpy.array (optional, default=None)
        weights given to ssm values. Required for weighted swi calculation.

    "next" refers to t_(i+1)
    "previous" refers to t_(i)

    Returns
    -------
    next_swi : numpy.array
        swi results
    next_den : numpy.array
        denominator for next iteration
    next_n : numpy.array
        total number of iterations next iteration
    next_wsum : numpy.array
        sum of weights for next iteration
    next_gain : numpy.array
        gain for next iteration
    """

    # Calculate time diff
    time_diff = (next_ssm_jd - jd).astype(np.float64)

    if weighted:

        next_e_term = np.exp(-(time_diff / float(tvalue)))

        # calculate denominator based on previous values
        next_den = (1.0 + next_e_term * den).astype(np.float32)

        # increase total number of iterations (i=0, regularly)
        next_n = (n + 1.0).astype(np.int32)

        # calculate next sum of weights
        next_wsum = (next_w + wsum).astype(np.float32)

        # calculate weighted SWI
        next_swi = ((swi * next_n / n * wsum * (next_den - 1) +
                     next_ssm * next_n * next_w) /
                    (next_wsum * next_den)).astype(np.float32)

        # rectify overshooting due to weighting
        next_swi[np.where(next_swi < 0)] = 0
        next_swi[np.where(next_swi > 100)] = 100

        return next_swi, next_den, next_n, next_wsum

    else:

        # calculate gain
        next_gain = (gain / (gain.astype(np.float64) +
                             np.exp(-(time_diff / float(tvalue)))))
        # calculate SWI
        next_swi = swi + next_gain * (next_ssm - swi)

        return next_swi, next_gain


class IterativeMultiSWI(object):
    """Class for calculation of the SWI in iterative mode.

    This means that previous_swi, previous_gain, previous_qflag and
    previous_julian_date have to be stored for a successful
    restart after a break.

    This class depends on the IterativeSWI class, but takes multiple T-values
    as argument.

    This is suitable if the SWI is calculated from a series of swaths/images.

    Parameters
    ----------
    ssm: numpy.ndarray
        An array of the soil moisture values.
    iter_data_path: string
        Path where the iteration data is/will be stored.
    tvalues: list
        List of tvalue integer values used for SWI calculation.
    weighted: boolean (default=False)
        Weighted swi calculation following equation (9) in the paper:
        Bauer-Marschallinger B, Paulik C, Hochstoeger S, Mistelbauer T,
        Modanesi S, Ciabatta L, Massari C, Brocca L, Wagner W.
        "Soil moisture from fusion of scatterometer and SAR:
        Closing the scale gap with temporal filtering."
        Remote Sensing. 2018 Jul;10(7):1030
    """

    def __init__(self, ssm, iter_data_path, tvalues, weighted=False):
        """
        Initializes the IterativeSWI class.

        Parameters
        ----------
        ssm: numpy.ndarray
            Soil Moisture values of the whole image.
        iter_data_path: string
            Output and input path of the calculation results.
        tvalues: list
            List of the tvalue integer values used for SWI calculation.
        weighted: boolean (default=False)
            Weighted swi calculation following equation (9) in the paper:
            Bauer-Marschallinger B, Paulik C, Hochstoeger S, Mistelbauer T,
            Modanesi S, Ciabatta L, Massari C, Brocca L, Wagner W.
            "Soil moisture from fusion of scatterometer and SAR:
            Closing the scale gap with temporal filtering."
            Remote Sensing. 2018 Jul;10(7):1030.
        """
        iter_swi = []
        prev_iter_files = []

        for i in range(0, len(tvalues)):
            it_swi = IterativeSWI(ssm, iter_data_path, tvalues[i], weighted=weighted)
            iter_swi.append(it_swi)
            prev_iter_files.append(it_swi.iterstepdata.files)

        self.calc_swi = iter_swi
        self.prev_iter_files = prev_iter_files
        self.tvalues = tvalues

    def load_iter_data(self):
        """
        Loads the iterative swi_img data for all tvalue values
        from the iter_data_path to iter_data.
        If there are no feasible files, the iter_data is set empty.
        """
        for it_data in self.calc_swi:
            it_data.load_iter_data()

    def store_iter_data(self):
        """ Stores the iterative swi_img data to the iter_data_path. """
        for it_data in self.calc_swi:
            it_data.store_iter_data()

    def remove_prev_iter_files(self):
        """ Removes the previous iter files in the iter_data_path. """
        for it_file_list in self.prev_iter_files:
            for it_file in it_file_list:
                os.remove(it_file)

    def calc_iter(self, next_ssm_jd, next_ssm, ind=None, weighted=False, next_w=None):
        """
        Calculate SWI iteratively.

        Parameters
        ----------
        next_ssm_jd: numpy.ndarray
            observation time of ssm given as julian date
        next_ssm: numpy.ndarray
            Soil moisture values
        ind: numpy.ndarray (optional, default=None)
            Index of valid next_ssm pixels.
            By default the complete image array is used.
        weighted : boolean (default=False)
            Weighted swi calculation following Equation (9) in the paper:
            Bauer-Marschallinger B, Paulik C, Hochstoeger S, Mistelbauer T,
            Modanesi S, Ciabatta L, Massari C, Brocca L, Wagner W.
            "Soil moisture from fusion of scatterometer and SAR:
            Closing the scale gap with temporal filtering."
            Remote Sensing. 2018 Jul;10(7):1030
        next_w : numpy.ndarray (optional, default=None)
            Weights given to ssm values.

        Returns
        ------
        swi_dict: dict
            dict with the new SWI values for all tvalue values.
            The keys of the dictionary are saved like e.g. SWI_001.
        """

        swi_dict = {}

        for i in range(0, len(self.tvalues)):
            key_swi = "SWI_%03d" % (self.tvalues[i],)
            swi_dict[key_swi] = self.calc_swi[i].calc_iter(next_ssm_jd, next_ssm, ind=ind,
                                                                weighted=weighted, next_w=next_w)

        return swi_dict


class IterativeSWI(object):
    """
    Class for calculation of the SWI in iterative mode.

    This means that the previous: swi, julian dates, qflags and/or gain, denominator, weights
    have to be stored for a successful restart after a break. On the first iteration the swi
    values are initialized as equal to the surface soil moisture, and gain is initialized as 1.0
    following Albergel et al (2008) https://doi.org/10.5194/hess-12-1323-2008

    This class only works with a single tvalue value.

    This is suitable if the SWI is calculated from a series of swaths/images.

    Parameters
    ----------
    ssm: numpy.ndarray
        Soil Moisture values of the whole image.
    iter_data_path: string
        Path to where the iteration data is/will be stored.
    tvalue: int
        T-value for SWI calculation
    weighted : boolean (default=False)
        Weighted swi calculation following equation (9) in the paper:
        Bauer-Marschallinger B, Paulik C, Hochstoeger S, Mistelbauer T,
        Modanesi S, Ciabatta L, Massari C, Brocca L, Wagner W.
        "Soil moisture from fusion of scatterometer and SAR:
        Closing the scale gap with temporal filtering."
        Remote Sensing. 2018 Jul;10(7):1030.
    """

    def __init__(self, ssm, iter_data_path, tvalue, weighted=False):
        """
        Initializes the IterativeSWI class.

        Parameters
        ----------
        ssm: numpy.ndarray
            Soil Moisture values of the whole image.
        iter_data_path: string
            Output and input path of the calculation results.
        tvalue: int
            T-value for SWI calculation
        weighted : boolean (default=False)
            Weighted swi calculation following equation (9) in the paper:
            Bauer-Marschallinger B, Paulik C, Hochstoeger S, Mistelbauer T,
            Modanesi S, Ciabatta L, Massari C, Brocca L, Wagner W.
            "Soil moisture from fusion of scatterometer and SAR:
            Closing the scale gap with temporal filtering."
            Remote Sensing. 2018 Jul;10(7):1030.
        """

        pref = "SWI_%03d" % (tvalue,)
        self.iter_data_path = iter_data_path
        self.float_nan = np.float32(np.nan)
        self.jd_nan = np.float64(np.nan)
        if weighted:
            self.iterstepdata = IterStepData(self.iter_data_path,
                                             len(ssm),
                                             {'swi': self.float_nan,
                                              'qflag': 1.0,
                                              'den': self.float_nan,
                                              'n': 0,
                                              'wsum': self.float_nan,
                                              'jd': self.jd_nan},
                                             prefix=pref)
        else:
            self.iterstepdata = IterStepData(self.iter_data_path,
                                             len(ssm),
                                             {'swi': self.float_nan,
                                              'qflag': 1.0,
                                              'gain': self.float_nan,
                                              'jd': self.jd_nan},
                                             prefix=pref)
        self.tvalue = tvalue
        self.processing_start = datetime.now()
        self.load_iter_data()

    def load_iter_data(self):
        """
        Loads the iterative swi_img data from the iter_data_path to iter_data.
        If there is no feasible file, the iter_data is set empty.
        """
        if self.iterstepdata.files_available is False:
            empty_data = self.iterstepdata.get_empty_data()
            empty_data['header']['processing_start'] = datetime(1900, 1, 1)
            empty_data['header']['processing_end'] = datetime(1901, 1, 1)
            self.iter_data = empty_data
        else:
            self.iter_data = self.iterstepdata.read_latest_iter_data()

    def store_iter_data(self):
        """
        Stores the iterative swi_img data to the iter_data_path.
        """
        self.iterstepdata.save_iter_data(self.iter_data)


    def calc_iter(self, next_ssm_jd, next_ssm, ind=None, weighted=False, next_w=None):
        """
        Calculate SWI iteratively.

        Parameters
        ----------
        next_ssm_jd: numpy.ndarray
            Observation time of ssm given as julian date
        next_ssm: numpy.ndarray
            Surface soil moisture values
        ind: numpy.ndarray
            Index of valid next_ssm pixels
        weighted: boolean (default=False)
            Weighted swi calculation following equation (9) in the paper:
            Bauer-Marschallinger B, Paulik C, Hochstoeger S, Mistelbauer T,
            Modanesi S, Ciabatta L, Massari C, Brocca L, Wagner W.
            "Soil moisture from fusion of scatterometer and SAR:
            Closing the scale gap with temporal filtering."
            Remote Sensing. 2018 Jul;10(7):1030.
        next_w : numpy.ndarray (optional, default=None)
            Weights given to ssm values

        Returns
        ------
        swi_img: numpy.ndarray
            Array with the new calculated swi_img values.
        """

        if ind is None:
            ind = np.arange(0, self.iter_data['swi'].size, 1)

        # load values from previous caluclation
        swi = self.iter_data['swi'][ind]
        jd = self.iter_data['jd'][ind]
        if weighted:
            den = self.iter_data['den'][ind]
            n = self.iter_data['n'][ind]
            wsum = self.iter_data['wsum'][ind]
        else:
            gain = self.iter_data['gain'][ind]

        still_untouched = np.isnan(swi)
        swi[still_untouched] = next_ssm[still_untouched]
        jd[still_untouched] = next_ssm_jd[still_untouched]
        if weighted:
            wsum[still_untouched] = next_w[still_untouched]
            den[still_untouched] = 1.0
            n[still_untouched] = 1
        else:
            gain[still_untouched] = 1.0

        # calculate new values
        if weighted:
            next_swi, next_den, next_n, next_wsum = calc_swi(next_ssm, next_ssm_jd,
                                       self.tvalue, swi, jd, den, n, wsum, next_w)
        else:
            next_swi, next_gain = calc_swi(next_ssm, next_ssm_jd, self.tvalue, swi, jd, gain=gain)

        # update values of iter data with now calculated values
        self.iter_data['swi'][ind] = next_swi
        self.iter_data['jd'][ind] = next_ssm_jd
        if weighted:
            self.iter_data['den'][ind] = next_den
            self.iter_data['n'][ind] = next_n
            self.iter_data['wsum'][ind] = next_wsum
        else:
            self.iter_data['gain'][ind] = next_gain

        # update iter_data header for complete file.
        valid_jd = np.isfinite(self.iter_data['jd'])
        header = {'sensing_start': julian2datetime(float(np.nanmin(self.iter_data['jd'][valid_jd]))),
                  'sensing_end': julian2datetime(float(np.max(self.iter_data['jd'][valid_jd]))),
                  'processing_start': self.processing_start,
                  'processing_end': datetime.now()}
        self.iter_data['header'].update(header)

        return self.iter_data['swi']
