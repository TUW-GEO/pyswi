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
from iterative_storage.iter_data import IterStepData
from pyswi.swi.swi import iterative_swi
from pytesmo.timedate.julian import julian2datetime
import numpy as np


class IterativeSWI(object):
    """Class for calculation of the SWI in iterative mode.

    This means that previous_swi, previous_gain, previous_qflag and
    previous_julian_date have to be stored for a successful
    restart after a break.

    This is suitable if the SWI is calculated from a series of swaths/images.

    Parameters
    ----------
    sm: numpy.recarray
        record array of the soil moisture values.
    iter_data_path: string
        Path where the iteration data is/will be stored.
    """

    def __init__(self, sm, iter_data_path):
        """
        Initializes the IterativeSWI class.

        Parameters
        ----------
        sm: numpy.ndarray
            Soil Moisture values of the whole image.
        iter_data_path: string
            Output and input path of the calculation results.

        """
        self.iter_data_path = iter_data_path
        self.float_nan = -999999.
        self.gain_nan = 1.
        self.swi_nan = 0.
        self.iterstepdata = IterStepData(self.iter_data_path,
                                         len(sm),
                                         {'swi': self.swi_nan,
                                          'gain': self.gain_nan,
                                          'qflag': self.float_nan,
                                          'jd': self.float_nan},
                                         prefix="SWI")
        self.processing_start = datetime.now()
        self.load_iter_data()

    def load_iter_data(self):
        """
        Loads the iterative swi data from the iter_data_path to iter_data.
        If there is no feasible file, the iter_data is set empty.
        """
        if self.iterstepdata.files_available is False:
            self.iter_data = self.iterstepdata.get_empty_data()
        else:
            self.iter_data = self.iterstepdata.read_latest_iter_data()

    def store_iter_data(self):
        """ Stores the iterative swi data to the iter_data_path. """
        self.iterstepdata.save_iter_data(self.iter_data)

    def calc_iter(self, sm_jd, sm, ctime):
        """
        Calculate SWI iteratively.

        Parameters
        ----------
        sm_jd: numpy.ndarray
            Julian days array
        sm: numpy.ndarray
            Soil moisture values
        ctime: int
            T value for SWI calculation

        Returns
        ------
        swi: numpy.ndarray
            Array with the new calculated swi values.
        """
        prev_swi = self.iter_data['swi']
        prev_gain = self.iter_data['gain']
        prev_qflag = self.iter_data['qflag']
        prev_jd = self.iter_data['jd']

        swi, qflag, gain = iterative_swi(sm, sm_jd, ctime, prev_swi, prev_gain,
                                         prev_qflag, prev_jd)

        self.iter_data['swi'] = swi
        self.iter_data['qflag'] = qflag
        self.iter_data['gain'] = gain
        self.iter_data['jd'] = sm_jd

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
