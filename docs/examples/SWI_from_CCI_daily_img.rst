This examples shows how to calculate SWI from daily ESA CCI SSM images using pyswi package.

Imports and dependencies. You will need to install the esa_cci_sm package first:
::

    import os
    import xarray as xr
    import numpy as np
    import pandas as pd
    from datetime import datetime
    from esa_cci_sm.interface import CCI_SM_025Img
    from pyswi.swi_img.iterative_swi import IterativeMultiSWI

List the desired T-values for the calculation, i.e.:
::

    t-values = [1, 5, 10, 15, 20, 40, 60, 100]

Read the first CCI SSM daily image (1d array for later indexing) and use it to initialize the SWI processor, path to store the results and a list of T-values for the SWI calculation.
::

    reader = CCI_SM_025Img(os.path.join(CCI_path, CCI_file))
    CCI_image = reader.read()
    SSM = CCI_image['sm'].flatten()

    swi_iter = IterativeMultiSWI(SSM, outpath, t_values)

Iterate over the CCI images, subseting each image to only valid SSM observations:
::

    ind = np.where(np.isfinite(SSM))[0]
    vSSM = SSM[ind]

Get the image's timestamp from datetime and convert it to julian date, then make it an array to match vSSM:
::

    dt = xr.open_dataset(os.path.join(CCI_path, CCI_file))
    dt = xr.DataArray(dt['time'].values)
    ts = pd.Timestamp(dt.values[0].astype(datetime))
    julian = ts.to_julian_date()
    jd = np.zeros_like(vSSM, dtype=np.float64)
    jd.fill(julian)

Calculate SWI for the supplied T-values:
::

    results = swi_iter.calc_iter(jd, vSSM, ind)