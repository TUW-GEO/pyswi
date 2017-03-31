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

from math import exp

uint8_nan = np.iinfo(np.uint8).max

float64_nan = -999999.
float32_nan = -999999.


def process_swi_pd(ssm_userformat_data, proc_param={},
                   ctime=[1, 5, 10, 15, 20, 40, 60, 100], gain_in=None,
                   jd_daily_out=False):
    """
    Processing of surface soil water index and the gain value.
    This function calls the calculate function and handles
    the conversions of the input/output formats of the gain to
    the formats that are needed for the calculation.

    Parameters
    ----------
    ssm_userformat_data : pandas.DataFrame
        Surface soil moisture userdata time series including ssm and ssf values.
    proc_param : dict
        Processing parameters. e.g. date_begin, date_end (julian date)
    ctime : list
        Ctime values that should be calculated.
    gain_in : pandas.DataFrame
         Input gain parameters of last calculation.
         fields gpi, last_jd, nom, denom
    jd_daily_out : bool
         If the flag is True the jd value of the swi will be set to every day
         23:59 of the time range, otherwise it takes the jd values from the ssm
         input within the range.

    Returns
    -------
    swi_ts : pandas.DataFrame
        Soil Water Index time series. Each field is for a ctime value.
        fields jd, SWI_xxx, QFLAG_xxx (e.g. SWI_010 for ctime value 10)
    gain : pandas.DataFrame
        Output gain parameters of last calculation.
        Each field is for a ctime value.
        fields gpi, last_jd(= date index), NOM_xxx, DENOM_xxx
        (e.g. NOM_010 for ctime value 10)
    """

    juldate = pd_datetime_to_juldate(ssm_userformat_data.axes[0].values)

    juldate = np.array(juldate)

    if gain_in is not None:
        gain_in = gain_dataframe_to_inputdict(gain_in, ctime)

    swi_ts, gain_out_mp = process_swi(ssm_userformat_data['sm'].values, juldate,
                                      proc_param=proc_param, ctime=ctime,
                                      gain_in=gain_in,
                                      jd_daily_out=jd_daily_out)

    swi_ts = swi_outputdict_to_dataframe(swi_ts, ctime)

    gain = gain_outputdict_to_dataframe(gain_out_mp, ctime)

    return swi_ts, gain


def process_swi(ssm, jd, proc_param={}, ctime=[1, 5, 10, 15, 20, 40, 60, 100],
                gain_in=None, jd_daily_out=False):
    """
    Processing of surface soil water index and the gain value.
    The ssm should already be filtered.
    This function calls the calculate function and handles
    the conversions of the input/output formats of the gain to
    the formats that are needed for the calculation.

    Parameters
    ----------
    ssm : numpy.ndarray
        Surface soil moisture time series.
    jd : numpy.ndarray
        Julian dates time series.
    proc_param : dict
        Processing parameters. e.g. date_begin, date_end (julian date)
    ctime : list
        Ctime values that should be calculated.
    gain_in : dict
         Input gain parameters of last calculation.
         fields gpi, last_jd, nom, denom
    jd_daily_out : bool
         If the flag is True the jd value of the swi will be set to every day
         23:59 of the time range, otherwise it takes the jd values from the ssm
         input within the range.

    Returns
    -------
    swi_ts : dict
        Soil Water Index time series. Each field is for a ctime value.
        fields jd, SWI_xxx, QFLAG_xxx (e.g. SWI_010 for ctime value 10)
    gain : dict
        Output gain parameters of last calculation.
        Each field is for a ctime value.
        fields gpi, last_jd, NOM_xxx, DENOM_xxx
        (e.g. NOM_010 for ctime value 10)
    """

    ssm_ts = {'jd': jd,
              'sm': ssm}

    swi_dates = []

    if 'date_from' in proc_param and 'date_to' in proc_param:
        date_from = proc_param['date_from']
        date_to = proc_param['date_to']

        if jd_daily_out is False:
            for d in range(0, len(jd)):
                if date_from <= jd[d] <= date_to:
                    swi_dates.append(jd[d])
            swi_dates = np.array(swi_dates)
    else:
        date_from = np.min(ssm_ts['jd'])
        date_to = np.max(ssm_ts['jd'])
        if jd_daily_out is False:
            swi_dates = jd

    num_swi = (date_to - date_from)
    num_swi = num_swi.astype('timedelta64[D]')
    num_swi = num_swi.astype(int)+1

    if jd_daily_out is True:
        month, day, year = julian.caldat(date_from)
        start_jd = julian.julday(month, day,
                                 year, 23, 59)

        swi_dates = np.arange(num_swi) + start_jd

    if type(ctime) is not np.array:
        ctime = np.array(ctime)

    # Init empty swi_ts
    swi_ts = {'jd': np.zeros([num_swi], dtype=np.float64),
              'swi': np.zeros([num_swi, len(ctime)], dtype=np.float32),
              'qflag': np.zeros([num_swi, len(ctime)], dtype=np.float32)}

    if num_swi != len(jd):
        indices = [i for i, x in enumerate(ssm_ts['jd']) if x == date_from]
        startdate_idx = indices[0]
        ssm_ts['sm'] = ssm[startdate_idx:]
        ssm_ts['jd'] = jd[startdate_idx:]

    swi_ts['jd'] = swi_dates

    gain = None

    if gain_in is not None:
        gain = inputdict_to_gain_dict(gain_in, ctime)

    swi_ts, gain_out_mp = calc(ssm_ts, swi_ts,
                               ctime=ctime, gain=gain)

    swit_ts = swi_dict_to_outputdict(swi_ts, ctime)

    gain_out = gain_dict_to_outputdict(gain_out_mp, ctime)

    return swit_ts, gain_out


def calc(ssm_ts, swi_ts, gain=None, ctime=np.array([1, 5, 10, 15, 20, 40, 60,
                                                    100])):
    """
    Calculation of surface soil water index.

    Parameters
    ----------
    ssm_ts : dict
        Surface soil moisture time series.
        fields sm, jd as numpy arrays.
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

    if gain is None:
        gain = {'last_jd': np.double(0),
                'denom': np.zeros(len(ctime)),
                'nom': np.zeros(len(ctime))
                }
        gain['denom'].fill(1)
        gain['nom'].fill(1)

    juldate = ssm_ts['jd']
    ssm = ssm_ts['sm']

    ssm = np.asarray(ssm, dtype=np.float32)

    norm_factor = \
        np.sum(np.exp(-(np.tensordot(np.arange(1000, dtype=np.float32),
                                     (1 / ctime), axes=0))), axis=0)

    last_jd_var = gain['last_jd']

    if gain['denom'][0] == 1:
        last_jd_var = juldate[0]

    nom = gain['nom']
    denom = gain['denom']

    swi_ts['swi'] = np.zeros([len(swi_ts['jd']), len(ctime)], dtype=np.float32)
    swi_ts['qflag'] = np.zeros([len(swi_ts['jd']), len(ctime)],
                               dtype=np.float32)

    swi_ts['swi'], swi_qflag, nom, denom, last_jd_var = \
        swi_calc_cy(juldate, ssm, ctime, swi_ts['jd'], swi_ts['swi'],
                    swi_ts['qflag'], nom, denom, last_jd_var, norm_factor)

    gain_out = {'denom': denom, 'nom': nom, 'last_jd': last_jd_var}

    return swi_ts, gain_out


def calc_noise(ssm_noise, jd, ctime=[1, 5]):
    """
    Calculation of surface soil water index noise.

    Parameters
    ----------
    ssm_noise : numpy.ndarray
        Surface Soil Moisture Noise time series.
    jd : numpy.ndarray
        Array of the julian dates.
    ctime : numpy.ndarray
        Integer values for the ctime variations.

    Returns
    -------
    swi_noise_ts : dict
        Soil water index noise series.
        fields jd, nom_sn, den_sn, swi_noise
    """
    len_sm = len(jd)
    len_ctime = len(ctime)

# var{n_swi} calculation
    nom_db = np.zeros((len_sm, len_ctime))
    den_db = np.zeros((len_sm, len_ctime))
    n_swi = np.zeros((len_sm, len_ctime))

    for i in range(0, len_sm):
        nom = np.zeros(len_ctime)
        den = np.zeros(len_ctime)
        for c in range(0, len_ctime):
            for j in range(0, i+1):
                exp_term = exp(((-2)*(jd[i]-jd[j])) / ctime[c])
                nom[c] = nom[c] + ssm_noise[i] * exp_term
                den[c] = den[c] + exp(((-1)*(jd[i] - jd[j])) / ctime[c])
            nom_db[i][c] = nom[c]
            den_db[i][c] = den[c]
            n_swi[i][c] = (nom[c] / (den[c] ** 2))

    inputdict = {'nom_sn': nom_db,
                 'den_sn': den_db,
                 'swi_noise': n_swi,
                 'jd': jd}

    variables_ctimedep = ['nom_sn', 'den_sn', 'swi_noise']
    variables_other = ['jd']

    swi_noise_ts = dict_to_outputdict(inputdict, variables_ctimedep, ctime,
                                      variables_other=variables_other)

    return swi_noise_ts


def calc_noise_rec(ssm_noise, jd, ctime=[1, 5]):
    """
    Calculation of surface soil water index noise.
    Recurrence approach of calculation.

    Parameters
    ----------
    ssm_noise : numpy.ndarray
        Surface Soil Moisture Noise time series.
    jd : numpy.ndarray
        Array of the julian dates.
    ctime : numpy.ndarray
        Integer values for the ctime variations.

    Returns
    -------
    swi_noise_ts : dict
        Soil water index noise series.
        fields jd, nom_sn, den_sn, swi_noise
    """
    len_sm = len(jd)
    len_ctime = len(ctime)

    # var{n_swi} calculation
    nom_sn = np.zeros((len_sm, len_ctime))
    den = np.zeros((len_sm, len_ctime))
    n_swi = np.zeros((len_sm, len_ctime))

    # exp(0) = 1
    nom_sn[0] = ssm_noise[0]
    den[0] = 1

    for i in range(0, len_sm):
        for c in range(0, len_ctime):
            fn = exp(((-1) * (jd[i]-jd[i-1]))/ctime[c])
            den[i][c] = 1 + fn * den[i-1][c]
            den_sn = den[i][c] ** 2

            exp_term = exp(((-2)*(jd[i]-jd[i-1])) / ctime[c])
            nom_sn[i][c] = nom_sn[i-1][c] * exp_term + ssm_noise[i]

            n_swi[i][c] = nom_sn[i][c] / den_sn

    inputdict = {'nom_sn': nom_sn,
                 'den_sn': den,
                 'swi_noise': n_swi,
                 'jd': jd}

    variables_ctimedep = ['nom_sn', 'den_sn', 'swi_noise']
    variables_other = ['jd']

    outputdict = dict_to_outputdict(inputdict, variables_ctimedep, ctime,
                                    variables_other=variables_other)

    return outputdict

# --- Help functions ---


def pd_datetime_to_juldate(pd_datetimes):
    """
    Converts panda.DatetimeIndex list format to a julian dates list.

    Parameters
    ----------
    pd_datetimes : numpy.ndarray
        Pandas Datetime array

    Returns
    -------
    juldates : numpy.ndarray
        Julian dates.
    """
    date_buffer = pd.DatetimeIndex(pd_datetimes)

    return date_buffer.to_julian_date()


def swi_outputdict_to_dataframe(swi_dict, ctime):
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

    jd = swi_dict['jd']
    dataframe_dict = {'jd': jd}

    for t in ctime:
        dataframe_dict["SWI_%03d" % (t,)] = swi_dict["SWI_%03d" % (t,)]
        dataframe_dict["QFLAG_%03d" % (t,)] = swi_dict["QFLAG_%03d" % (t,)]

    swi_ts_df = pd.DataFrame(dataframe_dict, index=dates)

    return swi_ts_df


def dict_to_outputdict(inputdict, variables_ctimedep, ctime,
                       variables_other=[]):
    """
    Converts a processed dictionary format to a usable output dictionary
    with names e.g. SWI_010 for swi as keys.

    Parameters
    ----------
    inputdict : dict
        Dictionary with variables that have different values per ctime.
        e.g. {swi : [[1,2,3], [1,2,3]]] for 3 different ctime values.
    variables_ctimedep: list
        List of variable names (inputdict key values) that have for every ctime
        different values
    ctime : numpy.ndarray
        Integer values for the ctime variations.
    variables_other: list
        Optional variables that are independent of the ctime values and
        are just added to the output dictionary without changes.

    Returns
    -------
    output_dict : dict
        Dictionary with ctime related key names for all variables_ctimedep.
        Each field is for a ctime value.
        (e.g. SWI_010 for ctime value 10)
    """

    depvar_buffer = {}

    for var in variables_ctimedep:
        depvar_buffer[var] = inputdict[var]

    output_dict = {}

    for var in variables_other:
        output_dict[var] = inputdict[var]

    for name, value in depvar_buffer.iteritems():
        for i in range(0, len(ctime)):
            if np.isscalar(value[i]):
                output_dict[name.upper()+"_%03d" % (ctime[i],)] = value[i]
            else:
                output_dict[name.upper()+"_%03d" % (ctime[i],)] = \
                    np.array([item[i] for item in value])

    return output_dict


def swi_dict_to_outputdict(swi_dict, ctime):
    """
    Converts a swi dictionary format to a usable output dictionary
    with names like SWI_010 as keys.

    Parameters
    ----------
    swi_dict : dict
        Soil Water Index time series as a dictionary.
        fields gpi, jd, swi, qflag
    ctime : numpy.ndarray
        Integer values for the ctime variations.

    Returns
    -------
    swi_ts_dict : dict
        Soil Water Index time series. Each field is for a ctime value.
        fields jd, SWI_xxx, QFLAG_xxx (e.g. SWI_010 for ctime value 10)
    """

    variables_ctimedep = ['swi', 'qflag']
    variables_other = ['jd']

    return dict_to_outputdict(swi_dict, variables_ctimedep, ctime,
                              variables_other=variables_other)


def gain_dict_to_outputdict(gain_dict, ctime):
    """
    Converts a gain(_out) dictionary format to an usable output format
    with keyvalues for each ctime.

    Parameters
    ----------
    gain_dict : dict
        Gain parameters of last calculation.
        fields gpi, last_jd, nom, denom
    ctime : numpy.ndarray
        Integer values for the ctime variations.

    Returns
    -------
    gain_dict : dict
        Gain parameters of last calculation. Each field is for a ctime value.
        fields last_jd, NOM_xxx, DENOM_xxx (e.g. NOM_010 for ctime value 10)
    """
    variables_ctimedep = ['nom', 'denom']
    variables_other = ['last_jd']

    return dict_to_outputdict(gain_dict, variables_ctimedep, ctime,
                              variables_other=variables_other)


def inputdict_to_gain_dict(gain_in, ctime):
    """
    Converts an input gain dictionary format to
    an algorithmusable gain dict format.
    There the ctime index equals the gain index fot that ctime.

    Parameters
    ----------
    gain_in : dict
        Gain parameters of last calculation.
        fields gpi, last_jd, nom, denom
    ctime : numpy.ndarray
        Integer values for the ctime variations.

    Returns
    -------
    gain_dict : dict
        Gain parameters of last calculation.
        fields last_jd, nom, denom
    """
    gain_dict = {'last_jd': gain_in['last_jd'], 'nom': np.zeros(len(ctime)),
                 'denom': np.zeros(len(ctime))}

    for i in range(0, len(ctime)):
        gain_dict['nom'][i] = gain_in["NOM_%03d" % (ctime[i],)]
        gain_dict['denom'][i] = gain_in["DENOM_%03d" % (ctime[i],)]

    return gain_dict


def gain_outputdict_to_dataframe(gain_dict, ctime):
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
        fields last_jd(= date index), NOM_xxx, DENOM_xxx
        (e.g. NOM_010 for ctime value 10)
    """
    last_jd = pd.to_datetime([gain_dict['last_jd'],
                              gain_dict['last_jd']] -
                              pd.Timestamp(0).to_julian_date(), unit='D')

    dataframe_dict = {}

    for i in range(0, len(ctime)):
        dataframe_dict["NOM_%03d" % (ctime[i],)] = \
            gain_dict["NOM_%03d" % (ctime[i],)]
        dataframe_dict["DENOM_%03d" % (ctime[i],)] = \
            gain_dict["DENOM_%03d" % (ctime[i],)]

    gain_df = pd.DataFrame(dataframe_dict, index=last_jd)

    return gain_df


def gain_dataframe_to_inputdict(gain_dataframe, ctime):
    """
    Converts a gain(_in) pandas dataframe to a dictionary format.

    Parameters
    ----------
    gain_dataframe : pandas.DataFrame
        Gain parameters of last calculation. Each field is for a ctime value.
        fields gpi, last_jd(= date index), NOM_xxx, DENOM_xxx
        (e.g. NOM_010 for ctime value 10)
    ctime : numpy.ndarray
        Integer values for the ctime variations.

    Returns
    -------
    gain_dict : dict
        Gain parameters of last calculation.
        fields last_jd, nom, denom
    """

    gain_dict = {}
    last_jd = gain_dataframe.index[0]
    gain_dict['last_jd'] = julian.julday(last_jd.month, last_jd.day,
                                         last_jd.year, last_jd.hour,
                                         last_jd.minute, last_jd.second)

    for i in range(0, len(ctime)):
        gain_dict["NOM_%03d" % (ctime[i],)] = \
            gain_dataframe["NOM_%03d" % (ctime[i],)].values[0]
        gain_dict["DENOM_%03d" % (ctime[i],)] = \
            gain_dataframe["DENOM_%03d" % (ctime[i],)].values[0]

    return gain_dict
