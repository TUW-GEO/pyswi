This examples shows how to calculate SWI time series from ASCAT surface soil moisture time series using pyswi package.

Import all necessary dependencies:

.. code:: python

    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore') # some import warnings are expected and ignored
        # install the ascat package first https://github.com/TUW-GEO/ascat
        from ascat.read_native.cdr import AscatSsmCdr

    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from numpy.lib.recfunctions import unstructured_to_structured
    from pyswi.swi_ts.swi_ts import calc_swi_ts

Set up the Ascat reader:

.. code:: python

    ascat_data_folder = os.path.join(testdata_path, 'sat', 'ascat', 'netcdf', '55R22')
    ascat_grid_folder = os.path.join(testdata_path, 'sat', 'ascat', 'netcdf', 'grid')
    static_layers_folder = os.path.join(testdata_path, 'sat', 'h_saf', 'static_layer')

    # init the AscatSsmCdr reader with the paths
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore') # some warnings are expected and ignored

        ascat_reader = AscatSsmCdr(ascat_data_folder, ascat_grid_folder,
                                   grid_filename='TUW_WARP5_grid_info_2_1.nc',
                                   static_layer_path=static_layers_folder)

Read a point in Northern Italy:

.. code:: python

    # point at (11,45)
    ascat_ts = ascat_reader.read(2302069, mask_ssf=True, mask_frozen_prob=80, mask_snow_prob=20)
    ascat_ts.plot()

.. image:: /_static/images/swi_calculation/output_1_1.png

Calculate SWI for disclosed T-values from soil moisture time series using swi_ts():

.. code:: python

    # Drop NA measurements
    ascat_sm_ts = ascat_ts.data[['sm']].dropna()

    # Get julian dates of time series
    jd = ascat_sm_ts.index.to_julian_date().values

    # T-values for SWI calculation
    t_values = [10, 50]

    sm_ts = {}
    sm_ts['jd'] = jd
    sm_ts['sm'] = ascat_sm_ts['sm'].values

    dtype = np.dtype([('jd', np.float64), ('sm', np.float32)])
    sm_ts = unstructured_to_structured(
        np.hstack((jd[:, np.newaxis], sm_ts['sm'][:, np.newaxis])), dtype=dtype)

    swi_ts, gain_out = calc_swi_ts(sm_ts, jd, t_value=t_values)

    fig, ax = plt.subplots(1, 1, figsize=(15,5))
    ax.plot(sm_ts['sm'], alpha=0.4, marker='o', color='#00bfff', label='SSM')
    ax.plot(swi_ts['swi_10'], lw=2, label='SWI T=10')
    ax.plot(swi_ts['swi_50'], lw=2, label='SWI T=50')
    plt.legend(loc=2)

.. image:: /_static/images/swi_calculation/output_2_1.png