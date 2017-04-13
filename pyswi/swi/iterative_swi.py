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
    """Class for calculation of the SSF in iterative mode.

    This means that mean_sigma, mean_sigma_gain, previous_ssf and
    previous_julian_date have to be stored for a successful restart after a break.

    This is suitable if the SSF is calculated from a series of swaths/images.

    Parameters
    ----------
    ft_thresholds: numpy.recarray
        record array of the freeze thaw thresholds of the full grid(mostly WARP DGG).

        with the fields:
        ['sig_flevel']
        ['sig_tlevel']
        ['pt1_doy']
        ['pt2_doy']
        ['msig_summer']
        ['msig_winter']
        ['sig_stdev_frozen']
        ['sig_slop_minust']
        ['pice_flag']
    iter_data_path: string
        Path where the iteration data is/will be stored.
    td_max: float, optional
        time difference maximum that is allowed between observation so the same points.
        This is given in fractional days. Default is 4 days since this is on the safe side for
        data in low latitude regions.
    """

    def __init__(self, sm, iter_data_path):
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
        if self.iterstepdata.files_available is False:
            self.iter_data = self.iterstepdata.get_empty_data()
        else:
            self.iter_data = self.iterstepdata.read_latest_iter_data()

    def store_iter_data(self):
        self.iterstepdata.save_iter_data(self.iter_data)

    def calc_iter(self, sm_jd, sm, ctime):

        prev_swi = self.iter_data['swi']
        prev_gain = self.iter_data['gain']
        prev_qflag = self.iter_data['qflag']
        prev_jd = self.iter_data['jd']

        swi, qflag, gain = iterative_swi(sm, sm_jd, ctime, prev_swi, prev_gain, prev_qflag, prev_jd)

        self.iter_data['swi'] = np.zeros([len(swi)])
        self.iter_data['gain'] = np.ones([len(swi)])
        self.iter_data['qflag'] = np.zeros([len(swi)])
        self.iter_data['jd'] = np.array([2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333, 2456596.64583333])
      #  self.iter_data['jd'].fill()
        #self.iter_data['swi'] = swi
        #self.iter_data['qflag'] = qflag
        #self.iter_data['gain'] = gain
        #self.iter_data['jd'] = sm_jd
        # update iter_data header for complete file.
        valid_jd = np.where(self.iter_data['jd'] != self.float_nan)
        header = {'sensing_start': julian2datetime(np.min(self.iter_data['jd'][valid_jd])),
                  'sensing_end': julian2datetime(np.max(self.iter_data['jd'][valid_jd])),
                  'processing_start': self.processing_start,
                  'processing_end': datetime.now()}

        self.iter_data['header'].update(header)
        return swi
