import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
import pylab as mpl
from meta import *
import matplotlib.dates as mdates

from cycler import cycler
import matplotlib

font = {'family': 'Arial',
        'size': 7}

matplotlib.rc('font', **font)
plt.rcParams['lines.linewidth'] = 1.0
matplotlib.rcParams['pdf.fonttype'] = 42

# mpl.rcParams['font.sans-serif'] = ['Times New Roman']  # font: SimSun，FangSong，SimHei，Times New Roman
# plt.rcParams['font.size'] = 8  # unit: pt
# plt.rcParams['lines.linewidth'] = 1  # line width
matplotlib.rcParams['axes.unicode_minus'] = False  # 
plt.rcParams['lines.linestyle'] = '-'  # - -- -. :
plt.rcParams['lines.markersize'] = 3  # marker size
plt.rcParams['lines.marker'] = 'None'  # None . , o v ^ < > s + x D d
plt.rcParams['mathtext.default'] = 'regular'  # 
plt.rcParams['mathtext.fontset'] = 'stix'  # formula font
# mpl.rcParams.update(mpl.rcParamsDefault) # default

plt.rc('axes', prop_cycle=(
    cycler(color=['#00796b', '#ec407a', '#6b81f2', '#2b821d', '#c2653c', '#801124', '#cda819', '#005eaa', '#0098d9'])))
plt.rc('axes', prop_cycle=(cycler(
    color=['#2ec7c9', '#b6a2de', '#5ab1ef', '#ffb980', '#d87a80', '#8d98b3', '#e5cf0d', '#c05050', '#95706d',
           '#dc69aa'])))

cls = ['#fcce10', '#b6a2de', '#5ab1ef', '#ffb980', '#d87a80', '#8d98b3', '#edafda', '#cbb0e3', '#f58db2', '#22c3aa', ]
fx_N = {'RNH4_ALGAE_UPTK': 'NH3 uptake by algae',
        'RNO3_ALGAE_UPTK': 'NO3 uptake by algae',
        'NH4_FROM_ALGAE': 'Metabolized NH3 by algae',
        'DON_FROM_ALGAE': 'Metabolized DON by algae',
        'PON_FROM_ALGAE': 'Metabolized PON by algae',
        'DON_FROM_PON': 'Hydrolysis of PON',
        'NH4_FROM_DON': 'Mineralization of DON ',
        'NO3_FROM_NH4': 'Nitrification of NH3', }
fx_P = {'PO4_ALGAE_UPTK': 'PO4 uptake by Algae',
        'PO4_FROM_ALGAE': 'Metabolized PO4 by algae',
        'DOP_FROM_ALGAE': 'Metabolized DOP by algae',
        'POP_FROM_ALGAE': 'Metabolized POP by algae',
        'DOP_FROM_POP': 'Hydrolysis of POP',
        'PO4_FROM_DOP': 'Mineralization of DOP', }


def read_result(scn):
    data = pd.read_excel(f'{flux_output_path}\\flux_NP_S{scn:02d}.xlsx')
    data['N_Deposition'] = data['ALGAEN_SET1'] + data['PON_SET1']
    data['Benthic_N'] = data['RNH4_BEN1'] + data['RNO3_BEN1']
    data['Inflow_N'] = data['ALGAEN_PS'] + data['PON_PS'] + data['DON_PS'] + data['RNH4_PS'] + data['RNO3_PS']
    data['Outflow_N'] = data['ALGAEN_EXIT'] + data['PON_EXIT'] + data['DON_EXIT'] + data['RNH4_EXIT'] + data[
        'RNO3_EXIT']
    data['total_input'] = data['Benthic_N'] + data['ALGAEN_fix'] + data['N_atm'] + data['N_resus'] + data['Inflow_N']
    data['total_output'] = data['Outflow_N'] + data['ALGAEN_SET1'] + data['PON_SET1'] + data['NO3_denit']
    data['P_Deposition'] = data['ALGAEP_SET1'] + data['POP_SET1'] + data['PO4_SET1']
    data['Benthic_P'] = data['PO4_BEN1']
    data['Inflow_P'] = data['ALGAEP_PS'] + data['POP_PS'] + data['DOP_PS'] + data['PO4_PS']
    data['Outflow_P'] = data['ALGAEP_EXIT'] + data['POP_EXIT'] + data['DOP_EXIT'] + data['PO4_EXIT']
    data['total_input'] = data['Benthic_P'] + data['P_atm'] + data['P_resus'] + data['Inflow_P']
    data['total_output'] = data['Outflow_P'] + data['ALGAEP_SET1'] + data['POP_SET1']
    return data


scn = 0
# Boundary Fluxes of N
data = read_result(scn)
data.loc[:, 'date'] = data['JDAY'].apply(lambda x: pd.Timedelta(days=x - 1) +
                                                   pd.to_datetime(jday_start))
dat_pos = data.copy()
dat_neg = data.copy()
float_cols = data.dtypes[data.dtypes.apply(lambda x: x == np.dtype('float64'))].index
dat_pos.loc[:, float_cols] = dat_pos.loc[:, float_cols].applymap(lambda x: max(x, 0))
dat_neg.loc[:, float_cols] = dat_neg.loc[:, float_cols].applymap(lambda x: min(x, 0))

fig, axes = plt.subplots(2, 2, figsize=(19/2.54, 9/2.54), sharex=True, dpi=500)
axes[0, 0].stackplot(dat_pos['date'], dat_pos['Inflow_N'], dat_pos['Outflow_N'], dat_pos['Benthic_N'],
                     dat_pos['ALGAEN_fix'],
                     dat_pos['N_atm'], dat_pos['N_resus'],
                     dat_pos['ALGAEN_SET1'], dat_pos['PON_SET1'], dat_pos['NO3_denit'],
                     colors=[cls[1], cls[9], cls[5], cls[8], cls[2], cls[3], cls[0], cls[6], cls[4]], alpha=1,
                     labels=['N inflow', 'N outflow', 'Benthic N exchange', 'N fixation',
                             'Atmospheric N deposition', 'N resuspension',
                             'Algal N sedimentation', 'PON sedimentation', 'Denitrification'])
axes[0, 0].stackplot(dat_neg['date'], dat_neg['Inflow_N'], dat_neg['Outflow_N'], dat_neg['Benthic_N'],
                     dat_neg['ALGAEN_fix'],
                     dat_neg['N_atm'], dat_neg['N_resus'],
                     dat_neg['ALGAEN_SET1'], dat_neg['PON_SET1'], dat_neg['NO3_denit'],
                     colors=[cls[1], cls[9], cls[5], cls[8], cls[2], cls[3], cls[0], cls[6], cls[4]], alpha=1, )

axes[0, 0].legend(loc='center left', bbox_to_anchor=(0.95, 0.5), frameon=False)
# axes[0, 0].set_xlabel('Time')
axes[0, 0].set(xlabel=None)
axes[0, 0].set_ylabel('Boundary fluxes of N $(t\ d^{-1})$')
axes[0, 0].spines['top'].set_visible(False)
axes[0, 0].spines['right'].set_visible(False)
axes[0, 0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
axes[0, 0].xaxis.set_major_locator(mdates.MonthLocator(interval=2))
axes[0, 0].tick_params(axis='x', labelrotation=30)
axes[0, 0].tick_params(axis='x', labelrotation=30)
axes[0, 0].text(0.05, 0.9, '(a)', transform=axes[0, 0].transAxes)

axes[0, 1].stackplot(dat_pos['date'], dat_pos['Inflow_P'], dat_pos['Outflow_P'],
                     dat_pos['ALGAEP_SET1'], dat_pos['POP_SET1'], dat_pos['Benthic_P'], dat_pos['P_atm'],
                     dat_pos['P_resus'],
                     colors=[cls[1], cls[9], cls[5], cls[4], cls[8], cls[2], cls[3]], alpha=1,
                     labels=['P inflow', 'P outflow', 'Algal P sedimentation', 'POP sedimentation',
                             'Benthic P exchange',
                             'Atmospheric P deposition', 'P resuspension',
                             ])
axes[0, 1].stackplot(dat_neg['date'], dat_neg['Inflow_P'], dat_neg['Outflow_P'],
                     dat_neg['ALGAEP_SET1'], dat_neg['POP_SET1'], dat_neg['Benthic_P'], dat_neg['P_atm'],
                     dat_neg['P_resus'],
                     colors=[cls[1], cls[9], cls[5], cls[4], cls[8], cls[2], cls[3]], alpha=1)

axes[0, 1].legend(loc='center left', bbox_to_anchor=(0.95, 0.5), frameon=False)
# axes[0, 1].set_xlabel('Time')
axes[0, 1].set(xlabel=None)
axes[0, 1].set_ylabel('Boundary fluxes of P $(t\ d^{-1})$')
axes[0, 1].spines['top'].set_visible(False)
axes[0, 1].spines['right'].set_visible(False)
axes[0, 1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
axes[0, 1].xaxis.set_major_locator(mdates.MonthLocator(interval=2))
axes[0, 1].tick_params(axis='x', labelrotation=30)
axes[0, 1].tick_params(axis='x', labelrotation=30)
axes[0, 1].text(0.05, 0.9, '(b)', transform=axes[0, 1].transAxes)
axes[0, 1].set_ylim([-4, 10])

for i, fx in enumerate(fx_N):
    axes[1, 0].plot(data['date'], data[fx], alpha=0.9, label=fx_N[fx])
axes[1, 0].legend(loc='center left', bbox_to_anchor=(0.95, 0.5), frameon=False)
# axes[1, 0].set_xlabel('Time')
axes[1, 0].set(xlabel=None)
axes[1, 0].set_ylabel('Intrawater fluxes of N $(t\ d^{-1})$')
axes[1, 0].spines['top'].set_visible(False)
axes[1, 0].spines['right'].set_visible(False)
axes[1, 0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
axes[1, 0].xaxis.set_major_locator(mdates.MonthLocator(interval=2))
axes[1, 0].tick_params(axis='x', labelrotation=30)
axes[1, 0].tick_params(axis='x', labelrotation=30)
axes[1, 0].text(0.05, 0.9, '(c)', transform=axes[1, 0].transAxes)

for i, fx in enumerate(fx_P):
    axes[1, 1].plot(data['date'], data[fx], alpha=0.9, label=fx_P[fx])
axes[1, 1].legend(loc='center left', bbox_to_anchor=(0.95, 0.5), frameon=False)
# axes[1, 1].set_xlabel('Time')
axes[1, 1].set(xlabel=None)
axes[1, 1].set_ylabel('Intrawater fluxes of P $(t\ d^{-1})$')
axes[1, 1].spines['top'].set_visible(False)
axes[1, 1].spines['right'].set_visible(False)
axes[1, 1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
axes[1, 1].xaxis.set_major_locator(mdates.MonthLocator(interval=2))
axes[1, 1].tick_params(axis='x', labelrotation=30)
axes[1, 1].tick_params(axis='x', labelrotation=30)
axes[1, 1].text(0.05, 0.9, '(d)', transform=axes[1, 1].transAxes)

fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
fig.tight_layout()
fig.savefig(f'{flux_stack_plot_path}\\stack_combined_S{scn:02d}_s2.jpg', dpi=300)
fig.savefig(f'{flux_stack_plot_path}\\stack_combined_S{scn:02d}_s2.pdf', dpi=500,
            format='pdf')
plt.show()

dat_pos.loc[:, 'Tot_N'] = (dat_pos['Inflow_N'] + dat_pos['Outflow_N'] + dat_pos['Benthic_N'] +
                           dat_pos['ALGAEN_fix'] + dat_pos['N_atm'] + dat_pos['N_resus'] +
                           dat_pos['ALGAEN_SET1'] + dat_pos['PON_SET1'] + dat_pos['NO3_denit'])
dat_pos.loc[:, 'Tot_P'] = (dat_pos['Inflow_P'] + dat_pos['Outflow_P'] + dat_pos['ALGAEP_SET1'] +
                           dat_pos['POP_SET1'] + dat_pos['Benthic_P'] + dat_pos['P_atm'] +
                           dat_pos['P_resus'])

dat_neg.loc[:, 'Tot_N'] = (dat_neg['Inflow_N'] + dat_neg['Outflow_N'] + dat_neg['Benthic_N'] +
                           dat_neg['ALGAEN_fix'] + dat_neg['N_atm'] + dat_neg['N_resus'] +
                           dat_neg['ALGAEN_SET1'] + dat_neg['PON_SET1'] + dat_neg['NO3_denit'])
dat_neg.loc[:, 'Tot_P'] = (dat_neg['Inflow_P'] + dat_neg['Outflow_P'] + dat_neg['ALGAEP_SET1'] +
                           dat_neg['POP_SET1'] + dat_neg['Benthic_P'] + dat_neg['P_atm'] +
                           dat_neg['P_resus'])
