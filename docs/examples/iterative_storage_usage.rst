This example shows the process of creating, storing, reading and updating iterative_storage data object.

Usage
=====

.. code:: python

   # create new iteration data object
   data_coll = IterStepData(self.itersavepath, 100,
                            variables={'gain_sigma': -9999.,
                                       'ssf': 255},
                            prefix="Iter_data")

   # get the empty array if no previous data exists and fill it with values
   if data_coll.files_available == False:
       iter_data = data_coll.get_empty_data()
   else:
       iter_data = data_coll.read_latest_iter_data()
   iter_data['ssf'][50] = 38
   iter_data['gain_sigma'][50] = 48

   # update the header
   iter_data['header']['sensing_start'] = datetime(2007, 01, 01, 12, 15, 22)
   iter_data['header']['sensing_end'] = datetime(2007, 01, 02, 11, 59, 58)
   iter_data['header']['processing_end'] = datetime(2014, 01, 01, 15, 20, 43)
   iter_data['header']['processing_start'] = datetime(2014, 01, 01, 14, 20, 43)
   data_coll.save_iter_data(iter_data)

   # This will create a file with the name
   # Iter_data_20070101121522Z_20070102115958Z_20140101152043Z.nc
   # The filename is created from the sensing_start, sensing_end and processing_end
   # fields in the header.

   restored_data = data_coll.read_latest_iter_data()
   # This reads the newest iteration data by sensing_end time.
