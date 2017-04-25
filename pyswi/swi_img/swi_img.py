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


def iterative_swi(ssm, ssm_jd, ctime, prev_swi, prev_gain, prev_qflag, prev_jd):
    """
    Takes input data and previous data and calculates swi_img
    values for next dates. All arrays have to be the same shape
    so that the whole batch can be calculated at once
    initializing and masking has to be done outside of this basic
    function

    Parameters
    ----------
    ssm : numpy.array
        surface soil moisture
    ssm_jd : numpy.array
        observation time of ssm given as julian date
    ctime : int
        T value for SWI calculation
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
    qflag = 1.0 + np.exp(-(time_diff / float(ctime))) * prev_qflag

    gain = (prev_gain / (prev_gain.astype(np.float) +
                         np.exp(-(time_diff / float(ctime)))))

    # calculate SWI
    swi = prev_swi + gain * (ssm - prev_swi)

    return swi, qflag, gain
