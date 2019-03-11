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

"""
This module represents the WARP Processing step Soil Water Index (SWI).
"""

import numpy as np


def iterative_swi(ssm, ssm_jd, tvalue, prev_swi, prev_gain, prev_qflag,
                  prev_jd):
    """
    Takes input data and previous data and calculates swi_img
    values for the next date. All arrays have to be the same shape
    so that the whole batch can be calculated at once
    initializing and masking has to be done outside of this basic
    function

    Parameters
    ----------
    ssm : numpy.array
        surface soil moisture
    ssm_jd : numpy.array
        observation time of ssm given as julian date
    tvalue : int
        T-value for SWI calculation
    prev_swi : numpy.array
        values of SWI at previous iteration
    prev_gain : numpy.array
        values of gain from previous iteration
    prev_qflag : numpy.array
        values of qflag from previous iteration
    prev_jd : numpy.array
        values of julian date from previous iteration

    Returns
    -------
    swi_img : numpy.array
        swi_img results
    qflag : numpy.array
        quality flag results
    gain : numpy.array
        gain for next iteration
    """
    # Calculate time diff
    time_diff = (ssm_jd - prev_jd).astype(np.float)

    # calculate qflag and gain based on previous values
    qflag = 1.0 + np.exp(-(time_diff / float(tvalue))) * prev_qflag

    gain = (prev_gain / (prev_gain.astype(np.float) +
                         np.exp(-(time_diff / float(tvalue)))))

    # calculate SWI
    swi = prev_swi + gain * (ssm - prev_swi)

    return swi, qflag, gain


def iterative_weighted_swi(next_ssm, next_w, next_ssm_jd, tvalue,
                           swi, den, n, wsum, jd):
    """
    Takes input data and existing data and
    calculates weighted SWI values for the next date.
    All arrays have to be the same shape so that the whole batch
    can be calculated at once. Initializing and masking has to be
    done outside of this basic function.

    Following Equation (9) in the paper...

    Bauer-Marschallinger B, Paulik C, Hochstoeger S, Mistelbauer T,
    Modanesi S, Ciabatta L, Massari C, Brocca L, Wagner W.
    "Soil moisture from fusion of scatterometer and SAR:
    Closing the scale gap with temporal filtering."
    Remote Sensing. 2018 Jul;10(7):1030.

    "next" refers to t_(i+1)
    "previous" refers to t_(i)

    Parameters
    ----------
    next_ssm : numpy.array
        surface soil moisture
    next_w : numpy.array
        weights given to ssm values
    next_ssm_jd : numpy.array
        observation time of ssm given as julian date
    tvalue : int
        T-value for SWI calculation
    swi : numpy.array
        values of SWI at previous iteration
    den : numpy.array
        values of denominator from previous iteration
    n : numpy.array
        values of total number of previous iterations
    wsum : numpy.array
        values of sum of weights applied in previous iteration
    jd : numpy.array
        values of julian date from previous iteration

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
    """

    # Calculate the exponential term
    time_diff = (next_ssm_jd - jd)
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

    return next_swi, next_den, next_n, next_wsum


def iterative_weighted_swi_qflag(next_ssm, next_w, next_ssm_jd, tvalue,
                                 swi, den, n, wsum, qflag, jd):
    """
    As iterative_weighted_swi(), but also calculates the QFLAG for each time step.
    This is not appropriate for non-constant time steps,
    but works for regular iteration intervals (e.g. daily calculation).

    "next" refers to t_(i+1)
    "previous" refers to t_(i)

    Parameters
    ----------
    next_ssm : numpy.array
        surface soil moisture
    next_w : numpy.array
        weights given to ssm values
    next_ssm_jd : numpy.array
        observation time of ssm given as julian date
    tvalue : int
        T-value for SWI calculation
    swi : numpy.array
        values of SWI at previous iteration
    den : numpy.array
        values of denominator from previous iteration
    n : numpy.array
        values of total number of previous iterations
    wsum : numpy.array
        values of sum of weights applied in previous iteration
    qflag : numpy.array
        values of qflag from previous iteration
    jd : numpy.array
        values of julian date from previous iteration

    Returns
    -------
    next_swi : numpy.array
        swi results
    next_qflag : numpy.array
        quality flag results
    next_den : numpy.array
        denominator for next iteration
    next_n : numpy.array
        total number of iterations next iteration
    next_wsum : numpy.array
        sum of weights for next iteration
    """

    # Calculate the exponential term
    time_diff = (next_ssm_jd - jd)
    next_e_term = np.exp(-(time_diff / float(tvalue)))

    # calculate qflag based on previous values
    next_qflag = (1.0 + next_e_term * qflag).astype(np.float32)

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

    return next_swi, next_qflag, next_den, next_n, next_wsum