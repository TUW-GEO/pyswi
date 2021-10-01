from swi_ts import calc_swi_noise_rec, calc_swi_noise, calc_swi_ts
import pandas as pd
import numpy as np
from numpy.lib.recfunctions import unstructured_to_structured
from c3s_sm.interface import C3STs
from xarray import CFTimeIndex
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

# read ISMN station list
ismn_path = r"R:\Projects\G3P\06_workspace\T-calibration\station_list\ismn_station_list.csv"
ismn_stations = pd.read_csv(ismn_path)

# path to C3S soil moisure data
c3s_path = r"R:\Projects\C3S_312b\07_data\v202012_TCDR\063_images_to_ts\combined-daily"

gpi = 319537 # full ts
# gpi = 384336

c3s_reader = C3STs(c3s_path, remove_nans=True)
c3s = c3s_reader.read(gpi)
c3s.index = CFTimeIndex(c3s.index).to_datetimeindex()
c3s = c3s.loc['2020-01-01':'2020-12-31']
c3s_reader.close()

ssm_ts = {}
ssm_ts['sm'] = c3s['sm']
ssm_ts['sm_noise'] = c3s['sm_uncertainty']
ssm_ts['jd'] = c3s.index.to_julian_date().values.astype(np.float64)

dtype = np.dtype([('jd', np.float64), ('sm', np.float32), ('sm_noise', np.float32)])
ssm_ts = unstructured_to_structured(np.hstack((ssm_ts['jd'][:, np.newaxis], ssm_ts['sm'][:, np.newaxis],
                                               ssm_ts['sm_noise'][:, np.newaxis])), dtype=dtype)

# dtype = np.dtype([('jd', np.float64), ('sm', np.float32)])
# ssm_ts = unstructured_to_structured(np.hstack((ssm_ts['jd'][:, np.newaxis], ssm_ts['sm'][:, np.newaxis])), dtype=dtype)


# swi_ts, gain_out = calc_swi_ts(ssm_ts, ssm_ts['jd'], gain_in=None, t_value=[5])


swi_ts, gain_out = calc_swi_noise(ssm_ts, t_value=[5])# swi_ts, gain_out = calc_swi_noise_rec(ssm_ts, t)

results = pd.DataFrame(data=swi_ts, index=c3s.index)
results['sm'] = ssm_ts['sm']
results['sm_noise'] = ssm_ts['sm_noise']

print("Bazinga!")