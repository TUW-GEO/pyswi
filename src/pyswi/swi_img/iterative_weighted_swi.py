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

import os
from datetime import datetime

import numpy as np

from iterative_storage.iter_data import IterStepData
from pytesmo.timedate.julian import julian2datetime

from pyswi.swi_img.calc import iterative_weighted_swi
from pyswi.swi_img.calc import iterative_weighted_swi_qflag


class IterativeMultiWeiSWI(object):
    """
    Class for calculation of the weighted SWI in iterative mode.

    This means that previous_swi, previous_gain, previous_qflag and
    previous_julian_date have to be stored for a successful
    restart after a break.

    This class depends on the IterativeSWI class, but takes multiple T-values
    as argument.

    This is suitable if the SWI is calculated from a series of swaths/images.

    Parameters
    ----------
    sm: numpy.recarray
        record array of the soil moisture values.
    iter_data_path: string
        Path where the iteration data is/will be stored.
    tvalues: numpy.array
        Array of the tvalue integer values used for the SWI calculation.

    """

    def __init__(self, ssm, iter_data_path, tvalues):
        """
        Initializes the IterativeSWI class.

        Parameters
        ----------
        ssm: numpy.ndarray
            Surface Soil Moisture values of the whole image.
        iter_data_path: string
            Output and input path of the calculation results.
        tvalues: numpy.array
            Array of the tvalue integer values used for the SWI calculation.
        """
        iter_swi = []
        prev_iter_files = []

        for i in range(0, len(tvalues)):
            it_swi = IterativeWeiSWI(ssm, iter_data_path, tvalues[i])
            iter_swi.append(it_swi)
            prev_iter_files.append(it_swi.iterstepdata.files)

        self.iterative_swi = iter_swi
        self.prev_iter_files = prev_iter_files
        self.tvalues = tvalues

    def load_iter_data(self):
        """
        Loads the iterative swi_img data for all tvalue values
        from the iter_data_path to iter_data.
        If there are no feasible files, the iter_data is set empty.
        """
        for it_data in self.iterative_swi:
            it_data.load_iter_data()

    def store_iter_data(self):
        """
        Stores the iterative swi_img data to the iter_data_path.
        """
        for it_data in self.iterative_swi:
            it_data.store_iter_data()

    def remove_prev_iter_files(self):
        """
        Removes the previous iter files in the iter_data_path.
        """
        for it_file_list in self.prev_iter_files:
            for it_file in it_file_list:
                os.remove(it_file)

    def calc_iter(self, next_ssm_jd, next_ssm, next_w, ind, calculate_qflag=False):
        """
        Calculate SWI iteratively.

        Parameters
        ----------
        next_ssm_jd : numpy.array
            observation time of ssm given as julian date
        next_ssm : numpy.array
            surface soil moisture
        next_w : numpy.array
            weights given to ssm values
        ind : numpy.array
            indices in iter_arr to which values of next_ssm correspond
        calculate_qflag : bool, optional
            should the qflag be calculated (only True if iteration over
            constant timesteps (e.g. daily)

        Returns
        ------
        swi_dict: dict
            dict with the new SWI values for all tvalue values.
            The keys of the dictionary are saved like e.g. SWI_050, GAIN_050.

        """
        swi_dict = {}

        for i in range(0, len(self.tvalues)):
            key = "SWI_%03d" % (self.tvalues[i],)
            swi_dict[key] = self.iterative_swi[i].calc_iter(
                next_ssm_jd, next_ssm, next_w, ind, calculate_qflag=calculate_qflag)

        return swi_dict


class IterativeWeiSWI(object):
    """
    Class for calculation of the SWI in iterative mode.

    This means that previous_swi, previous_gain, previous_qflag and
    previous_julian_date have to be stored for a successful
    restart after a break.
    This class only works with one tvalue value.

    This is suitable if the SWI is calculated from a series of swaths/images.

    Parameters
    ----------
    sm: numpy.recarray
        record array of the soil moisture values.
    iter_data_path: string
        Path where the iteration data is/will be stored.
    tvalue: int
        T-value for the SWI calculation
    """

    def __init__(self, ssm, iter_data_path, tvalue):
        """
        Initializes the IterativeSWI class.

        Parameters
        ----------
        ssm: numpy.ndarray
            Soil Moisture values of the whole image.
        iter_data_path: string
            Output and input path of the calculation results.
        """
        pref = "SWI_%03d" % (tvalue,)
        self.iter_data_path = iter_data_path
        self.float_nan = np.float32(np.nan)
        self.iterstepdata = IterStepData(self.iter_data_path,
                                         len(ssm),
                                         {'swi': self.float_nan,
                                          'qflag': 1.0,
                                          'den': self.float_nan,
                                          'n': 0,
                                          'wsum': self.float_nan,
                                          'jd': np.float64(np.nan)},
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

    def calc_iter(self, next_ssm_jd, next_ssm, next_w, ind, calculate_qflag=False):
        """
        Calculate weighted SWI iteratively.

        Parameters
        ----------
        next_ssm_jd : numpy.array
            observation time of ssm given as julian date
        next_ssm : numpy.array
            surface soil moisture
        next_w : numpy.array
            weights given to ssm values
        ind : numpy.array
            indices in iter_arr to which values of next_ssm correspond
        calculate_qflag : bool, optional
            should the qflag be calculated (only True if iteration over
            constant timesteps (e.g. daily)

        Returns
        ------
        swi_img: numpy.ndarray
            Array with the new calculated swi_img values.
        """
        # load values from previous calculation
        swi = self.iter_data['swi'][ind]
        den = self.iter_data['den'][ind]
        n = self.iter_data['n'][ind]
        wsum = self.iter_data['wsum'][ind]
        jd = self.iter_data['jd'][ind]
        if calculate_qflag:
            qflag = self.iter_data['qflag'][ind]

        still_untouched = np.isnan(swi)
        swi[still_untouched] = next_ssm[still_untouched]
        wsum[still_untouched] = next_w[still_untouched]
        den[still_untouched] = 1.0
        n[still_untouched] = 1
        jd[still_untouched] = next_ssm_jd[still_untouched]
        if calculate_qflag:
            qflag[still_untouched] = 1.0

        # calculate new values
        if calculate_qflag:
            next_swi, next_qflag, next_den, next_n, next_wsum = \
                iterative_weighted_swi_qflag(next_ssm, next_w, next_ssm_jd,
                                       self.tvalue, swi, den, n, wsum, qflag, jd)
        else:
            next_swi, next_den, next_n, next_wsum = \
                iterative_weighted_swi(next_ssm, next_w, next_ssm_jd,
                                       self.tvalue, swi, den, n, wsum, jd)


        # update values of iter_data with now calculated values
        self.iter_data['swi'][ind] = next_swi
        self.iter_data['den'][ind] = next_den
        self.iter_data['n'][ind] = next_n
        self.iter_data['wsum'][ind] = next_wsum
        self.iter_data['jd'][ind] = next_ssm_jd
        if calculate_qflag:
            self.iter_data['qflag'][ind] = next_qflag

        # update iter_data header for complete file.
        valid_jd = np.isfinite(self.iter_data['jd'])
        header = {'sensing_start': julian2datetime(np.min(self.iter_data['jd'][valid_jd])),
                  'sensing_end': julian2datetime(np.max(self.iter_data['jd'][valid_jd])),
                  'processing_start': self.processing_start,
                  'processing_end': datetime.now()}

        self.iter_data['header'].update(header)
        return self.iter_data['swi']
