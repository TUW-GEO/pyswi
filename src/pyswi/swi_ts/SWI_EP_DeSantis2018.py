import numpy as np
from math import sqrt

def swi_ep_rec(ssm, T_value, T_noise):
    """
    Recursive SWI error propagation function based on DeSantis and Biondi (2018)
    "Error propagation from remotely sensed surface soil moisture
    into soil water index using an exponential filter",
    https://doi.org/10.29007/kvhb
    Translated from MatLab code obtained from the authors.

    Parameters
    ----------
    ssm: numpy.ndarray
        Surface soil moisture time series with fields 'sm', 'sm_uncertainty' and 'sm_jd'
    T_value: numpy.ndarray
        Exponential filter characteristic T-value parameter
    T_noise: numpy.ndarray
        T-value standard error. Baseed on empirical experiments
        author suggests 10% of T for calibrated T-values.

    Returns
    -------
    swi: numpy.ndarray
        Soil water index time series
    swi_noise:numpy.ndarray
        Soil water index noise time series
    """

    len_ssm = len(ssm)
    len_T = len(T_value)

    swi = np.zeros((len_ssm, len_T))
    swi_noise = np.zeros((len_ssm, len_T))

    swi[0] = ssm['sm'][0]
    swi_noise[0] = ssm['sm_uncertainty'][0]

    last_jd = ssm['sm_jd'][0]
    gain_curr = 1
    contr1_curr = ssm['sm_uncertainty'][0]**2
    G_curr = 0
    JT_curr = 0 # Jacobian term

    for i in range(1, len_ssm):
        for c in range(0, len_T):
            time_diff = ssm['sm_jd'][i]-last_jd
            Esponenz = np.exp(-time_diff / T_value[c]) # exp_term?
            gain_old = gain_curr
            contr1_old = contr1_curr
            gain_curr = gain_old / (gain_old + Esponenz)
            swi[i][c] = swi[i-1][c] + gain_curr * (ssm['sm'][i] - swi[i-1][c])
            contr1_curr = ((1 - gain_curr)**2) * contr1_old + (gain_curr * ssm['sm_uncertainty'][i])**2
            G_old = G_curr
            JT_old = JT_curr
            G_curr = Esponenz * (G_old + time_diff / T_value[c] / gain_old)
            JT_curr = gain_curr / T_value[c] * (G_curr * (swi[i-1][c] - swi[i][c])
                                                + Esponenz * T_value[c] / gain_old * JT_old)
            contr2 = (JT_curr * T_noise[c])**2
            swi_noise[i][c] = sqrt(contr1_curr + contr2)

        last_jd = ssm['sm_jd'][i]

    # todo:
    # gain_out = {} #does it make sense to restart time series?
    # it would mean reprocessing ts for a gpi, as in making an ICDR ts, and restarting this
    # this would always mean img -> ts first, so why not implement directly in IterativeSWI()

    dtype_list = [('swi_jd', np.float64)]
    for t in T_value:
        dtype_list.append(('swi_noise_{}'.format(t), np.float32))
        dtype_list.append(('swi_{}'.format(t), np.float32))

    swi_ts = np.zeros(ssm.size, dtype=np.dtype(dtype_list))

    swi_ts['swi_jd'] = ssm['sm_jd']
    for i, t in enumerate(T_value):
        swi_ts['swi_{}'.format(t)] = swi[:, i]
        swi_ts['swi_noise_{}'.format(t)] = swi_noise[:, i]

    return swi_ts #,gain_out for restarting? rather implement in iterative?


