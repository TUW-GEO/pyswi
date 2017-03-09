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
This module represents the WARP Processing step Soil Water Index (SWI).
"""
import numpy as np
import pyximport
pyximport.install(setup_args={'include_dirs': [np.get_include()]})
from pyswi.swi.swi_calc_routines import swi_calc_cy
import pytesmo.timedate.julian as julian
import pandas as pd

uint8_nan = np.iinfo(np.uint8).max

float64_nan = -999999.
float32_nan = -999999.


def process(ssm_userformat_data, proc_param, swi_param, gain_in=None):
    """
    Processing of surface soil water index and the gain value.
    This function calls the calculate function and handles the conversions of the input/output formats of the gain to
    the formats that are needed for the calculation.

    Parameters
    ----------
    ssm_userformat_data : pandas.DataFrame
        Surface soil moisture userdata time series including ssm and ssf values.
    proc_param : dict
        Processing parameters.
    swi_param : dict
        Soil Water Index parameters.
    gain_in : pandas.DataFrame
         Input gain parameters of last calculation.
         fields gpi, last_jd, nom, denom

    Returns
    -------
    swi_ts : pandas.DataFrame
        Soil Water Index time series. Each field is for a ctime value.
        fields jd, SWI_xxx, QFLAG_xxx (e.g. SWI_010 for ctime value 10)
    gain : pandas.DataFrame
        Output gain parameters of last calculation. Each field is for a ctime value.
        fields gpi, last_jd(= date index), NOM_xxx, DENOM_xxx (e.g. NOM_010 for ctime value 10)
    """

    ssm_ts = {'jd': ssm_userformat_data.axes[0].values,
              'sm': ssm_userformat_data['sm'].values}

    ssf_ts = {'ssf': ssm_userformat_data['ssf'].values}

    if 'data_from' in proc_param and 'date_to' in proc_param:
        date_from = proc_param['date_from']
        date_to = proc_param['date_to']
    else:
        date_from = np.min(ssm_ts['jd'])
        date_to = np.max(ssm_ts['jd'])

    num_swi = (date_to - date_from)
    num_swi = num_swi.astype('timedelta64[D]')
    num_swi = num_swi.astype(int)+1

    if swi_param['ctime'] is None:
        ctime = np.array([1, 5, 10, 15, 20, 40, 60, 100])
    else:   # Parsing list items from String
        ctime = swi_param['ctime']
        ctime = ctime[1:-1]
        ctime = np.fromstring(ctime, dtype=int, sep=',')

    # Init empty swi_ts
    swi_ts = {'jd': np.zeros([num_swi], dtype=np.float64),
              'swi': np.zeros([num_swi, len(ctime)], dtype=np.float32),
              'qflag': np.zeros([num_swi, len(ctime)], dtype=np.float32)}

    date_buffer = pd.DatetimeIndex([date_from])

    start_date = julian.julday(date_buffer.month[0], date_buffer.day[0], date_buffer.year[0], 00, 01)

    swi_dates = np.arange(num_swi) + start_date

    swi_ts['jd'] = swi_dates

    gain = None

    if gain_in is not None:
        gain = gain_dataframe_to_dict(gain_in, ctime)

    swi_ts, gain_out_mp = calc(ssm_ts, ssf_ts, swi_ts, ctime=ctime, gain=gain)

    swi_ts = swi_dict_to_dataframe(swi_ts, ctime)

    gain = gain_dict_to_dataframe(gain_out_mp, ctime)

    return swi_ts, gain


def calc(ssm_ts, ssf_ts, swi_ts, gain=None, ctime=None):
    """
    Calculation of surface soil water index.

    Parameters
    ----------
    ssm_ts : dict
        Surface soil moisture time series.
        fields sm, jd
    ssf_ts : dict
        Surface state flag time series.
        fields ssf
    swi_ts : dict
        Empty SWI time series dictionary.
    gain : dict
        Gain parameters of last calculation.
        fields gpi, last_jd, nom, denom
    ctime : numpy.ndarray
        Integer values for the ctime variations.

    Returns
    -------
    swi_ts : dict
        Soil water index series.
        fields jd, swi, qflag
    gain_out : dict
        Gain parameters of last calculation.
        fields gpi, last_jd, nom, denom
    """

    if ctime is None:
        ctime = np.array([1, 5, 10, 15, 20, 40, 60, 100])

    if gain is None:
        gain = {'last_jd': np.double(0),
                'denom': np.zeros(len(ctime)),
                'nom': np.zeros(len(ctime))
                }
        gain['denom'].fill(1)
        gain['nom'].fill(1)

#   Only considering "non frozen" values (ssf = 1)
    f_nan = float32_nan

    frozen = np.where(ssf_ts['ssf'] != 1)[0]

    ssm_ts['sm'][frozen] = f_nan

    valid = np.where(ssm_ts['sm'] != f_nan)[0]

    if len(valid) == 0:
        return swi_ts, gain

    juldate = ssm_ts['jd'][valid]
    ssm = ssm_ts['sm'][valid]
    ssm = np.asarray(ssm, dtype=np.float32)

    date_buffer = pd.DatetimeIndex(juldate)

    juldate_tmp = np.zeros([len(juldate)], dtype=np.float64)
    for i in range(0, len(date_buffer)):
        buf = julian.julday(date_buffer.month[i], date_buffer.day[i], date_buffer.year[i], date_buffer.hour[i],
                            date_buffer.minute[i], date_buffer.second[i])
        juldate_tmp[i] = buf

    juldate = juldate_tmp

    norm_factor = \
        np.sum(np.exp(-(np.tensordot(np.arange(1000, dtype=np.float32),
                                     (1 / ctime), axes=0))), axis=0)

#   last_jd_var is the index value of the gain file
    last_jd_var = gain['last_jd']

    if gain['denom'][0] == 1:
        last_jd_var = juldate[0]

    nom = gain['nom']
    denom = gain['denom']

    swi_ts['swi'] = np.zeros([len(swi_ts['jd']), len(ctime)], dtype=np.float32)
    swi_ts['qflag'] = np.zeros([len(swi_ts['jd']), len(ctime)], dtype=np.float32)
#   round juldate to make it possible to check equality.
    jd_round = swi_ts['jd'].round(2)
    juldate_round = juldate.round(2)
    last_jd_var = last_jd_var.round(2)

    swi_ts['swi'], swi_qflag, nom, denom, last_jd_var = swi_calc_cy(juldate_round, ssm, ctime, jd_round, swi_ts['swi'],
                                                                    swi_ts['qflag'], nom, denom, last_jd_var,
                                                                    norm_factor)

    gain_out = {'denom': denom, 'nom': nom, 'last_jd': last_jd_var}

    return swi_ts, gain_out

# --- Help functions ---


def swi_dict_to_dataframe(swi_dict, ctime):
    """
    Converts a swi dictionary format to a pandas DataFrame format.

    Parameters
    ----------
    swi_dict : dict
        Soil Water Index time series as a dictionary.
        fields gpi, jd, swi, qflag
    ctime : numpy.ndarray
        Integer values for the ctime variations.

    Returns
    -------
    swi_ts_df : pandas.DataFrame
        Soil Water Index time series. Each field is for a ctime value.
        fields jd, SWI_xxx, QFLAG_xxx (e.g. SWI_010 for ctime value 10)
    """
    dates = pd.to_datetime(swi_dict['jd'] - pd.Timestamp(0).to_julian_date(),
                           unit='D')

    swi = swi_dict['swi']
    jd = swi_dict['jd']
    qflag = swi_dict['qflag']

    dataframe_dict = {'jd': jd}

    for i in range(0, len(ctime)):
        dataframe_dict["SWI_%03d" % (ctime[i],)] = [item[i] for item in swi]
        dataframe_dict["QFLAG_%03d" % (ctime[i],)] = [item[i] for item in qflag]

    swi_ts_df = pd.DataFrame(dataframe_dict, index=dates)

    return swi_ts_df


def gain_dict_to_dataframe(gain_dict, ctime):
    """
    Converts a gain(_out) dictionary format to a pandas DataFrame format.

    Parameters
    ----------
    gain_dict : dict
        Gain parameters of last calculation.
        fields gpi, last_jd, nom, denom
    ctime : numpy.ndarray
        Integer values for the ctime variations.

    Returns
    -------
    gain_dataframe : pandas.DataFrame
        Gain parameters of last calculation. Each field is for a ctime value.
        fields gpi, last_jd(= date index), NOM_xxx, DENOM_xxx (e.g. NOM_010 for ctime value 10)
    """
    last_jd = pd.to_datetime([gain_dict['last_jd'], gain_dict['last_jd']] - pd.Timestamp(0).to_julian_date(), unit='D')

    nom = gain_dict['nom']
    denom = gain_dict['denom']

    dataframe_dict = {}

    for i in range(0, len(ctime)):
        dataframe_dict["NOM_%03d" % (ctime[i],)] = nom[i]
        dataframe_dict["DENOM_%03d" % (ctime[i],)] = denom[i]

    gain_df = pd.DataFrame(dataframe_dict, index=last_jd)

    return gain_df


def gain_dataframe_to_dict(gain_dataframe, ctime):
    """
    Converts a gain(_in) pandas dataframe to a dictionary format.

    Parameters
    ----------
    gain_dataframe : pandas.DataFrame
        Gain parameters of last calculation. Each field is for a ctime value.
        fields gpi, last_jd(= date index), NOM_xxx, DENOM_xxx (e.g. NOM_010 for ctime value 10)
    ctime : numpy.ndarray
        Integer values for the ctime variations.

    Returns
    -------
    gain_dict : dict
        Gain parameters of last calculation.
        fields gpi, last_jd, nom, denom
    """

    gain_dict = {}
    last_jd = gain_dataframe.index[0]
    gain_dict['last_jd'] = julian.julday(last_jd.month, last_jd.day, last_jd.year, last_jd.hour,
                                         last_jd.minute, last_jd.second)
    gain_dict['nom'] = np.zeros(len(ctime))
    gain_dict['denom'] = np.zeros(len(ctime))

    for i in range(0, len(ctime)):
        gain_dict['nom'][i] = gain_dataframe["NOM_%03d" % (ctime[i],)].values[0]
        gain_dict['denom'][i] = gain_dataframe["DENOM_%03d" % (ctime[i],)].values[0]

    return gain_dict
