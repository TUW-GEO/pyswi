# Copyright (c) 2021, Vienna University of Technology (TU Wien),
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
Storage of Iteration data between processing runs.
'''
import numpy as np
import os
import glob
from datetime import datetime
import netCDF4
from trollsift import Parser


class IterStepData(object):

    """
    reads and saves data needed between iteration steps

    The two methods get_empty_data and save_iter_data can be overwritten
    if more complex data structures are necessary.
    By default one dimensional arrays of size n_obs and a simple header are used.

    The default header is:

        header = {'product_name': "",
                  'parent_product_name': "",
                  'instrument_id': "",
                  'sensing_start': datetime(1900, 1, 1),
                  'sensing_end': datetime(4100, 12, 31, 15, 0, 0),
                  'processing_end': ""}

    Parameters
    ----------
    path : string
        path to which the iteration step data should be read from
        and saved to
    n_obs: int
        number of observations to store for the full target grid
    variables: dict
        a mapping from variable name to nan/fill value to use.
        This sets the datatype and initial values for the variable
        that is stored.
    prefix: string, optional
        filename prefix
    """

    def __init__(self, path, n_obs,
                 variables,
                 prefix='GENERIC_PREFIX'):
        self.path = path
        self.n_obs = n_obs
        self.variables = variables
        self.prefix = prefix
        self.files = 0
        self.metadata = None
        self.datestoragefmt = "%Y-%m-%d %H:%M:%S"
        self.fname_datefmt = '%Y%m%d%H%M%SZ'
        self.fname_template = '_'.join(
            [self.prefix,
             "{{sensing_start:{}}}".format(self.fname_datefmt),
             "{{sensing_end:{}}}".format(self.fname_datefmt),
             "{{processing_end:{}}}.nc".format(self.fname_datefmt)])
        self.fname_parser = Parser(self.fname_template)
        self._load_info()

    def _construct_file_name(self, sensing_start, sensing_end,
                             processing_end):
        """
        construct filename from sensing start and end and
        processing end

        Parameters
        ----------
        sensing_start : datetime
            time of observation of first measurement in file
        sensing_end : datetime
            time of observation of last measurement in file
        processing_end : datetime
            processing timestamp

        Returns
        -------
        filename : string
            filename
        """

        return self.fname_parser.compose({'sensing_start': sensing_start,
                                          'sensing_end': sensing_end,
                                          'processing_end': processing_end})

    def _load_info(self):
        """
        load the information about the files from the path
        """
        self.files = self.find_files()
        if len(self.files) == 0:
            self.files_available = False
        else:
            self.files_available = True
        self.metadata = self.get_metadata()

    def find_files(self):
        """
        find files in path by searching for 'self.prefix*'

        Returns
        -------
        files : list
            list of full filepaths
        """
        files = glob.glob(os.path.join(self.path, ''.join([self.prefix, '*'])))
        files.sort()

        return files

    def get_metadata(self):
        """
        loads metadata from filenames. Important field is the end of sensing
        in each file which indicates if a bufr file has already been processed
        """
        metadata = {'sensing_start': [],
                    'sensing_end': [],
                    'processing_end': []}

        fmt = '%Y%m%d%H%M%S'

        for _file in self.files:
            filename = os.path.split(_file)[1]
            mdata = self.fname_parser.parse(filename)
            metadata['sensing_start'].append(mdata['sensing_start'])
            metadata['sensing_end'].append(mdata['sensing_end'])
            metadata['processing_end'].append(mdata['processing_end'])

        return metadata

    def _get_latest_idx(self):
        """
        get index of latest sensing file
        """
        if type(self.metadata) != dict:
            self.get_metadata()

        max_sensing_idx = np.argmax(self.metadata['sensing_end'])

        return max_sensing_idx

    def get_latest_sensing_time(self):
        """
        find the latest time of the observations the iteration files in the
        folder already contain
        """

        if type(self.metadata) != dict:
            self.get_metadata()

        max_sensing_time = self.metadata['sensing_end'][self._get_latest_idx()]

        return max_sensing_time

    def read_latest_iter_data(self):
        """
        reads the latest iteration data
        """

        filename = self.files[self._get_latest_idx()]
        structure = {'header': {}}
        with netCDF4.Dataset(filename, 'r') as nc:
            for variable in nc.variables:
                try:
                    structure[variable] = nc.variables[variable][:].filled()
                except AttributeError:
                    structure[variable] = nc.variables[variable][:]
            for name in nc.ncattrs():
                structure['header'][name] = getattr(nc, name)

            from_string = ['sensing_start', 'sensing_end',
                           'processing_start', 'processing_end']
            for name in from_string:
                structure['header'][name] = datetime.strptime(
                    structure['header'][name], self.datestoragefmt)

        return structure

    def get_empty_data(self):
        """
        returns a empty structure dictionary
        """

        header = {'product_name': "",
                  'parent_product_name': "",
                  'instrument_id': "",
                  'sensing_start': datetime(1900, 1, 1),
                  'sensing_end': datetime(4100, 12, 31, 15, 0, 0),
                  'processing_end': ""}

        structure = {'header': header}
        for variable in self.variables:
            nan_array = np.empty(
                self.n_obs, dtype=type(self.variables[variable]))
            nan_array.fill(self.variables[variable])
            structure[variable] = nan_array

        return structure

    def save_iter_data(self, data, refresh_coll=True):
        """
        save iteration data to file and refresh the collection if set

        Parameters
        ----------
        data : dict
            dictionary in the format given by get_empty_data
        refresh_coll : boolean
            if True then the collection is refreshed after the file is saved
            which is good if processing should continue with new iteration
            data
        """
        filename = self._construct_file_name(data['header']['sensing_start'],
                                             data['header']['sensing_end'],
                                             data['header']['processing_end'])

        with netCDF4.Dataset(os.path.join(self.path, filename), 'w') as nc:
            nc.createDimension('observations', self.n_obs)
            header = data['header'].copy()
            # convert date object to string for storage in netCDF
            to_string = ['sensing_start', 'sensing_end',
                         'processing_start', 'processing_end']
            for name in to_string:
                header[name] = header[name].strftime(self.datestoragefmt)

            nc.setncatts(header)
            for variable in self.variables:
                var = nc.createVariable(variable, data[variable].dtype, ('observations'),
                                        fill_value=self.variables[variable], zlib=True)
                var[:] = data[variable]

        if refresh_coll:
            self._load_info()
