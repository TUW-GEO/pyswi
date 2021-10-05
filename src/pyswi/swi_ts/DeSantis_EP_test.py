from swi_ts import calc_swi_noise_rec, calc_swi_noise, calc_swi_ts
import pandas as pd
import numpy as np
from numpy.lib.recfunctions import unstructured_to_structured
from c3s_sm.interface import C3STs
from xarray import CFTimeIndex
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

from SWI_EP_DeSantis2018 import swi_ep_rec

# path to C3S soil moisure data
c3s_path = r"R:\Projects\C3S_312b\07_data\v202012_TCDR\063_images_to_ts\combined-daily"

gpi = 319537 # full ts
# gpi = 384336

c3s_reader = C3STs(c3s_path, remove_nans=True)
c3s = c3s_reader.read(gpi)
c3s.index = CFTimeIndex(c3s.index).to_datetimeindex()
c3s = c3s.loc['2020-01-01':'2020-12-31']
c3s_reader.close()

# ssm = {}
# ssm['sm'] = c3s['sm']
# ssm['sm_uncertainty'] = c3s['sm_uncertainty']
# ssm['sm_jd'] = c3s.index.to_julian_date().values.astype(np.float64)
#
# dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32), ('sm_uncertainty', np.float32)])
# ssm = unstructured_to_structured(np.hstack((ssm['sm_jd'][:, np.newaxis], ssm['sm'][:, np.newaxis],
#                                                ssm['sm_uncertainty'][:, np.newaxis])), dtype=dtype)
#
# T_values = [5]
# T_noise = [.5]
#
# swi = swi_ep_rec(ssm, T_values, T_noise)
#
# results = pd.DataFrame(data=swi, index=c3s.index)
# results['sm_uncertainty'] = ssm['sm_uncertainty']
# results['sm'] = ssm['sm']

ssm = {}
ssm['sm'] = c3s['sm']
ssm['sm_noise'] = c3s['sm_uncertainty']
ssm['sm_uncertainty'] = c3s['sm_uncertainty']
ssm['jd'] = c3s.index.to_julian_date().values.astype(np.float64)

dtype = np.dtype([('jd', np.float64), ('sm', np.float32), ('sm_noise', np.float32)])
ssm = unstructured_to_structured(np.hstack((ssm['jd'][:, np.newaxis], ssm['sm'][:, np.newaxis],
                                               ssm['sm_noise'][:, np.newaxis])), dtype=dtype)

swi_old, gain_out = calc_swi_ts(ssm, ssm['jd'], t_value=[5, 10])

ssm = {}
ssm['sm'] = c3s['sm']
ssm['sm_uncertainty'] = c3s['sm_uncertainty']
ssm['sm_jd'] = c3s.index.to_julian_date().values.astype(np.float64)

dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32), ('sm_uncertainty', np.float32)])
ssm = unstructured_to_structured(np.hstack((ssm['sm_jd'][:, np.newaxis], ssm['sm'][:, np.newaxis],
                                               ssm['sm_uncertainty'][:, np.newaxis])), dtype=dtype)

swi_new = swi_ep_rec(ssm, T_value=[5, 10], T_noise=[.5, 1.])

print("Bazinga!")

plt.plot(ssm['sm'])
plt.plot(swi_old['swi_5'])
plt.plot(swi_new['swi_5'])
plt.legend(['sm', 'swi_old_5', 'swi_new_5'], fontsize=20, loc='best')

#
# plt.plot(results['sm_uncertainty'])
# plt.plot(results['swi_noise_5'])
# plt.plot(results['calc_swi_noise'])
# plt.plot(results['calc_swi_noise_rec'])
# plt.plot(results['swi_ts_noise'], c='black')
#
# plt.legend(['sm_uncertainty', 'DeSantis&Biondi 2018 (\u03C3T=10%)', 'calc_swi_noise()',
# 'calc_swi_noise_rec()', 'calc_swi_ts()'], fontsize=20, loc='best')
#
# plt.title("Comparison of SWI noise caluclating functions, T=5", size=20)
#
# plt.ylabel('uncertainty [m3/m3]', fontsize=20, labelpad=20)
# plt.yticks(fontsize=20)
# plt.ylim((0.001, 0.008))
# plt.xticks(fontsize=20)


fig, [ax1, ax2] = plt.subplots(1,2)

ax1.plot(results['swi_noise_5'])

ax2.plot()

fig.legend(['DeSantis&Biondi 2018 (\u03C3T=10%)'])
print("Bazinga!")