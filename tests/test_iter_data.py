# Copyright (c) 2017, Vienna University of Technology (TU Wien),
# Department of Geodesy and Geoinformation (GEO).
# All rights reserved.
#
# All information contained herein is, and remains the property of Vienna
# University of Technology (TU Wien), Department of Geodesy and Geoinformation
# (GEO). The intellectual and technical concepts contained herein are
# proprietary to Vienna University of Technology (TU Wien), Department of
# Geodesy and Geoinformation (GEO). Dissemination of this information or
# reproduction of this material is forbidden unless prior written permission
# is obtained from Vienna University of Technology (TU Wien), Department of
# Geodesy and Geoinformation (GEO).

'''
Test for the iteration data storage
'''

import os
import unittest
import glob
import tempfile
import numpy.testing as nptest
import numpy as np
from datetime import datetime

from pyswi.swi_img.iterative_storage import IterStepData


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth


class TestIterStepData(unittest.TestCase):

    def setUp(self):

        self.itersavepath = tempfile.mkdtemp()

    def test_metadata_extraction(self):
        path = os.path.join(curpath(), 'data', 'iterative_storage', 'metadata_from_filename')
        data_coll = IterStepData(path, 100, {})

        assert data_coll.metadata['sensing_start'] == [
            datetime(2013, 11, 24, 12)]
        assert data_coll.metadata['sensing_end'] == [
            datetime(2013, 11, 24, 12, 2, 58)]
        assert data_coll.metadata['processing_end'] == [
            datetime(2013, 11, 25, 17, 47, 33)]

    def test_get_latest_sensing_time(self):
        path = os.path.join(curpath(), 'data', 'iterative_storage', 'metadata_from_filename')
        data_coll = IterStepData(path, 100, {})

        max_sensing_time = data_coll.get_latest_sensing_time()
        assert max_sensing_time == datetime(2013, 11, 24, 12, 2, 58)

    def test_construct_filename(self):
        path = os.path.join(curpath(), 'data', 'iterative_storage', 'metadata_from_filename')
        data_coll = IterStepData(path, 100, {}, prefix='ASCA_SWI25_03_M02')

        filename = data_coll._construct_file_name(datetime(2013, 11, 24, 12),
                                                  datetime(
                                                      2013, 11, 24, 12, 2, 58),
                                                  datetime(2013, 11, 25, 17, 47, 33))

        assert filename == 'ASCA_SWI25_03_M02_20131124120000Z_20131124120258Z_20131125174733Z.nc'

    def test_iter_data_read_write(self):

        data_coll = IterStepData(self.itersavepath, 100,
                                 variables={'gain_sigma': -9999.,
                                            'ssf': 255})
        iter_data = data_coll.get_empty_data()
        assert iter_data['ssf'][0] == 255
        assert iter_data['ssf'].dtype == int
        assert iter_data['gain_sigma'][0] == -9999.
        assert iter_data['gain_sigma'].dtype == np.float64
        iter_data['ssf'][50] = 38
        iter_data['gain_sigma'][50] = 48
        iter_data['header']['sensing_start'] = datetime(2007, 1, 1, 12, 15, 22)
        iter_data['header']['sensing_end'] = datetime(2007, 1, 2, 11, 59, 58)
        iter_data['header']['processing_end'] = datetime(2014, 1, 1, 15, 20, 43)
        iter_data['header']['processing_start'] = datetime(2014, 1, 1, 14, 20, 43)
        data_coll.save_iter_data(iter_data)
        restored_data = data_coll.read_latest_iter_data()
        restored_header = restored_data.pop('header')
        iter_header = iter_data.pop('header')
        for name in iter_header:
            assert restored_header[name] == iter_header[name]

        for name in restored_header:
            assert restored_header[name] == iter_header[name]

        for name in iter_data:
            nptest.assert_allclose(restored_data[name], iter_data[name])
