from src.pyswi.swi_ts.swi_ts import calc_swi_noise_rec, calc_swi_ts, swi_error_prop
import pandas as pd
import numpy as np
from numpy.lib.recfunctions import unstructured_to_structured
from c3s_sm.interface import C3STs
from xarray import CFTimeIndex
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

####################
t_value = [5, 50, 100]

swi_jd = pd.date_range(
    '2007-01-01', periods=9).to_julian_date().values.astype(np.float64)

sm = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90], dtype=np.float32)

dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32)])
ssm_ts = unstructured_to_structured(
    np.hstack((swi_jd[:, np.newaxis], sm[:, np.newaxis])), dtype=dtype)

swi_ts, gain_out = calc_swi_ts(ssm_ts, swi_jd, t_value=t_value)

swi_ref = {
    'swi_5': np.array([10., 15.49834, 21.32452, 27.472094, 33.93228, 40.69421,
                       47.7452, 55.07107, 62.65647]),
    'swi_50': np.array([10., 15.049998, 20.133324, 25.249971, 30.399931, 35.58319,
                        40.799732, 46.049545, 51.332603]),
    'swi_100': np.array([10., 15.025, 20.066666, 25.124996, 30.199991, 35.29165,
                         40.399967, 45.524944, 50.666576])
}

for t in t_value:
    np.testing.assert_array_almost_equal(swi_ts['swi_{}'.format(t)],
                                         swi_ref['swi_{}'.format(t)], 4)
##############

# path to C3S soil moisure data
c3s_path = r"R:\Projects\C3S_312b\07_data\v202012_TCDR\063_images_to_ts\combined-daily"

gpi = 319537 # full ts
# gpi = 384336

c3s_reader = C3STs(c3s_path, remove_nans=True)
c3s = c3s_reader.read(gpi)
c3s.index = CFTimeIndex(c3s.index).to_datetimeindex()
c3s = c3s.loc['2020-01-01':'2020-12-31']
c3s1 = c3s.loc['2020-01-01':'2020-06-30']
# c3s1 = c3s.loc['2020-01-01':'2020-06-29']
c3s2 = c3s.loc['2020-07-01':'2020-12-31']
c3s_reader.close()

c3s1['sm_uncertainty']['2020-06-30'] = np.nan

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
ssm['sm_jd'] = c3s.index.to_julian_date().values.astype(np.float64)

dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32), ('sm_noise', np.float32)])
ssm = unstructured_to_structured(np.hstack((ssm['sm_jd'][:, np.newaxis], ssm['sm'][:, np.newaxis],
                   ssm['sm_noise'][:, np.newaxis])), dtype=dtype)

ssm1 = {}
ssm1['sm'] = c3s1['sm']
ssm1['sm_noise'] = c3s1['sm_uncertainty']
ssm1['sm_jd'] = c3s1.index.to_julian_date().values.astype(np.float64)

dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32), ('sm_noise', np.float32)])
ssm1 = unstructured_to_structured(np.hstack((ssm1['sm_jd'][:, np.newaxis], ssm1['sm'][:, np.newaxis],
                   ssm1['sm_noise'][:, np.newaxis])), dtype=dtype)

ssm2 = {}
ssm2['sm'] = c3s2['sm']
ssm2['sm_noise'] = c3s2['sm_uncertainty']
ssm2['sm_jd'] = c3s2.index.to_julian_date().values.astype(np.float64)

dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32), ('sm_noise', np.float32)])
ssm2 = unstructured_to_structured(np.hstack((ssm2['sm_jd'][:, np.newaxis], ssm2['sm'][:, np.newaxis],
                   ssm2['sm_noise'][:, np.newaxis])), dtype=dtype)

# swi_old, gain_out = calc_swi_ts(ssm, ssm['jd'], t_value=[10])

# ssm = {}
# ssm['sm'] = c3s['sm']
#
# ssm['sm_uncertainty'] = c3s['sm_uncertainty']
# ssm['sm_jd'] = c3s.index.to_julian_date().values.astype(np.float64)
#
# dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32), ('sm_uncertainty', np.float32)])
# ssm = unstructured_to_structured(np.hstack((ssm['sm_jd'][:, np.newaxis], ssm['sm'][:, np.newaxis],
#                                                ssm['sm_uncertainty'][:, np.newaxis])), dtype=dtype)
######################## swi_error_prop
# swi_new = swi_error_prop(ssm, T_value=[5], T_noise=[.5])
# swi_multi, gain_out = swi_error_prop(ssm, t_value=[5,10], t_noise=[.5, 1])
# df = pd.DataFrame(swi_multi, columns=['swi_jd', 'swi_noise_5', 'swi_5', 'swi_noise_10', 'swi_10'], index=c3s.index)
#
# swi_split1, gain_out1 = swi_error_prop(ssm1, t_value=[5, 10], t_noise=[.5, 1])
# df1 = pd.DataFrame(swi_split1, columns=['swi_jd', 'swi_noise_5', 'swi_5', 'swi_noise_10', 'swi_10'], index=c3s1.index)
#
# swi_split2, gain_out2 = swi_error_prop(ssm2, gain_in=gain_out1, t_value=[5, 10], t_noise=[.5, 1])
# df2 = pd.DataFrame(swi_split2, columns=['swi_jd', 'swi_noise_5', 'swi_5', 'swi_noise_10', 'swi_10'], index=c3s2.index)
#
# swi_split3, gain_out3 = swi_error_prop(ssm2, t_value=[5, 10], t_noise=[.5, 1])
# df3 = pd.DataFrame(swi_split3, columns=['swi_jd', 'swi_noise_5', 'swi_5', 'swi_noise_10', 'swi_10'], index=c3s2.index)
#

####################### calc_swi_ts
swi_multi, gain_out = calc_swi_ts(ssm, ssm['sm_jd'], t_value=[5,10])
df = pd.DataFrame(swi_multi, columns=['swi_jd', 'swi_noise_5', 'swi_5', 'swi_noise_10', 'swi_10'], index=c3s.index)

swi_split1, gain_out1 = calc_swi_ts(ssm1, ssm1['sm_jd'], gain_in=None, t_value=[5,10])
df1 = pd.DataFrame(swi_split1, columns=['swi_jd', 'swi_noise_5', 'swi_5', 'swi_noise_10', 'swi_10'], index=c3s1.index)

swi_split2, gain_out2 = calc_swi_ts(ssm2, ssm2['sm_jd'], gain_in=gain_out1, t_value=[5,10])
df2 = pd.DataFrame(swi_split2, columns=['swi_jd', 'swi_noise_5', 'swi_5', 'swi_noise_10', 'swi_10'], index=c3s2.index)

print("Bazinga!")

df['swi_noise_5'].plot()
df1['swi_noise_5'].plot()
df2['swi_noise_5'].plot()

# plt.plot(ssm['sm'])
plt.plot(ssm['sm_uncertainty'])
# plt.plot(swi_new['swi_noise_5'])
plt.legend(['C3S SSM uncertainty', 'DeSantis&Biondi2018 (t=5, \u03C3T=10%)'], fontsize=20, loc='best')

plt.title("SWI (T=5) error propagation for a point in Southeastern Australia", fontsize=20)
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
plt.ylabel('uncertainty [m3/m3]', fontsize=20, labelpad=20)
plt.xlabel('day of year [2020]', fontsize=20, labelpad=20)
# plt.ylim((0.001, 0.008))

plt.yticks(fontsize=20)
plt.xticks(fontsize=20)


fig, [ax1, ax2] = plt.subplots(1,2)

ax1.plot(results['swi_noise_5'])

ax2.plot()

fig.legend(['DeSantis&Biondi 2018 (\u03C3T=10%)'])
print("Bazinga!")