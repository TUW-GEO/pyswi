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

from datetime import datetime

import numpy as np
from iterative_storage.iter_data import IterStepData
from pytesmo.timedate.julian import julian2datetime

from pyswi.swi_img.calc import iterative_weighted_swi


class IterativeMultiWeiSWI(object):
    """Class for calculation of the weighted SWI in iterative mode.

    This means that previous_swi, previous_gain, previous_qflag and
    previous_julian_date have to be stored for a successful
    restart after a break.

    This class depends on the IterativeSWI class, but takes multiple tvalues as
    argument.

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

        for i in range(0, len(tvalues)):
            it_swi = IterativeWeiSWI(ssm, iter_data_path, tvalues[i])
            iter_swi.append(it_swi)

        self.iterative_swi = iter_swi
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
        """ Stores the iterative swi_img data to the iter_data_path. """
        for it_data in self.iterative_swi:
            it_data.store_iter_data()

    def calc_iter(self, next_ssm_jd, next_ssm, next_w):
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

        Returns
        ------
        swi_dict: dict
            dict with the new SWI values for all tvalue values.
            The keys of the dictionary are saved like e.g. SWI_050, GAIN_050.

        """
        swi_dict = {}

        for i in range(0, len(self.tvalues)):
            key = "SWI_%03d" % (self.tvalues[i],)
            swi_dict[key] = self.iterative_swi[i].calc_iter(next_ssm_jd, next_ssm, next_w)

        return swi_dict


class IterativeWeiSWI(object):
    """Class for calculation of the SWI in iterative mode.

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
        self.float_nan = -999999.
        self.iterstepdata = IterStepData(self.iter_data_path,
                                         len(ssm),
                                         {'swi': self.float_nan,
                                          'qflag': self.float_nan,
                                          'den': self.float_nan,
                                          'n': self.float_nan,
                                          'wsum': self.float_nan,
                                          'jd': self.float_nan},
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
            self.iter_data = self.iterstepdata.get_empty_data()
        else:
            self.iter_data = self.iterstepdata.read_latest_iter_data()

    def store_iter_data(self):
        """ Stores the iterative swi_img data to the iter_data_path. """
        self.iterstepdata.save_iter_data(self.iter_data)

    def calc_iter(self, next_ssm_jd, next_ssm, next_w):
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

        Returns
        ------
        swi_img: numpy.ndarray
            Array with the new calculated swi_img values.

        """

        # load values from previous calculation
        swi = self.iter_data['swi']
        qflag = self.iter_data['qflag']
        den = self.iter_data['den']
        n = self.iter_data['n']
        wsum = self.iter_data['wsum']
        jd = self.iter_data['jd']

        # calculate new values
        next_swi, next_qflag, next_den, next_n, next_wsum = \
            iterative_weighted_swi(next_ssm, next_w, next_ssm_jd,
                                   self.tvalue, swi, den, n, wsum, qflag, jd)

        # update values of iter_data with now calculated values
        self.iter_data['swi'] = next_swi
        self.iter_data['qflag'] = next_qflag
        self.iter_data['den'] = next_den
        self.iter_data['n'] = next_n
        self.iter_data['wsum'] = next_wsum
        self.iter_data['jd'] = next_ssm_jd

        # update iter_data header for complete file.
        valid_jd = np.where(self.iter_data['jd'] != self.float_nan)
        header = {'sensing_start': julian2datetime(
                      np.min(self.iter_data['jd'][valid_jd])),
                  'sensing_end': julian2datetime(
                      np.max(self.iter_data['jd'][valid_jd])),
                  'processing_start': self.processing_start,
                  'processing_end': datetime.now()}

        self.iter_data['header'].update(header)
        return swi