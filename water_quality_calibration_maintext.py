import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pylab import mpl
import matplotlib

from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from math import sqrt
from math import log

from meta import *
from obs_compare import read_result

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


stations = {'Huiwanzhong': (11, 55), 'Luojiaying': (18, 49), 'GYSdong': (26, 40), 'GYSzhong': (19, 36),
            'GYSxi': (14, 34), 'Baiyukou': (16, 30), 'Haikouxi': (14, 19), 'Dianchinan': (21, 12)}


colors = ['k', 'r', 'g']
alpha = 0.8
colors = ['#194571', '#c51f27', 'g']
alpha = 1.0
bounds = {'TEM': 0.1, 'TN': 0.2, 'TP': 0.2, 'COD': 0.15, 'NH4': 0.2, 'CHLA': 0.3, 'DO': 0.15, 'ELE': 0.0001}
indexes = ['ELE', 'TEM', 'DO', 'TN', 'TP', 'CHLA']
units = {'ELE': 'Water level (m)', 'TEM': 'Water temperature ($^\mathrm{o}$C)', 'DO': 'DO (mg/L)',
         'TN': 'TN (mg/L)', 'TP': 'TP (mg/L)', 'CHLA': 'Chla (mg/L)'}

alldata = read_result('../rerun/S00/WQWCTS.OUT')
alldata = read_result('../baseline_DC-0107/WQWCTS.OUT')
alldata = read_result('D:/Files/Projects/NonLinearity/S00-12-18/WQWCTS.OUT')
observeds = pd.read_excel('.\data\observed.xlsx')
observeds.loc[:, 'DAY'] = observeds.DAY - 365
observeds = observeds.loc[(observeds.DAY > 0) & (observeds.DAY <= 2556)]

ELE = pd.read_excel('data/waterLevel.xlsx')
ELE.loc[:, 'DAY'] = ELE.DAY - 365
ELE = ELE.loc[(ELE.DAY > 0) & (ELE.DAY <= 2556)]

wq_simu, wq_obs = {}, {}
station = 'Huiwanzhong'
simulate = alldata.loc[
    (alldata['TIME'] % 1 > 0.45) & (alldata['TIME'] % 1 < 0.55) &
    (alldata['I'] == stations[station][0]) & (alldata['J'] == stations[station][1]),
    ['TIME', 'TEM', 'TN', 'TP', 'COD', 'NH4', 'CHLA', 'DO', 'ELE']]
observed = observeds.loc[(observeds['I'] == stations[station][0]) &
                         (observeds['J'] == stations[station][1])].copy()
simulate['DATE'] = pd.date_range('2012-1-1', '2018-12-30')
simulate.index = simulate['TIME'].astype(int)
observed.index = observed['DAY']
wq_simu[station] = simulate
wq_obs[station] = observed

fig, axes = plt.subplots(2, 3, figsize=(14 / 2.54, 7 / 2.54), sharex=True)
# plt.tick_params(labelsize=12)

for i, j, index in zip([0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2], indexes):

    axes[i, j].fill_between(wq_simu[station]['TIME'], wq_simu[station][index] * (1 - bounds[index]),
                            wq_simu[station][index] * (1 + bounds[index]), alpha=0.5, facecolor=colors[2])
    line1, = axes[i, j].plot(wq_simu[station]['TIME'], wq_simu[station][index], c=colors[0], alpha=alpha,
                             label='Simulated data')

    if index == 'ELE':
        line2, = axes[i, j].plot(ELE['DAY'], ELE[index], 'o', c=colors[1],
                                 markersize=2, markeredgewidth=0,
                                 alpha=alpha, label='Observed data')
        ymax = max((wq_simu[station][index] * (1 + bounds[index])).max(), ELE[index].max())
        ymax = 1890
        axes[i, j].set_ylim(1884, 1890)
    else:
        line2, = axes[i, j].plot(wq_obs[station]['DAY'], wq_obs[station][index], 'o', c=colors[1],
                                 markersize=2, markeredgewidth=0,
                                 alpha=alpha, label='Observed data')
        ymax = max((wq_simu[station][index] * (1 + bounds[index])).max(), wq_obs[station][index].max())
        if index == 'DO':
            axes[i, j].set_ylim(4, 14)
        else:
            axes[i, j].set_ylim(0, ymax)

    axes[i, j].vlines(x=0.6935, ymin=0, ymax=1, color='#EC5D3B', linestyles='dashed',
                      transform=axes[i, j].transAxes)
    # axes[i, j].set_title(station)
    # axes[i,j].set_yticks(list(np.linspace(0,ymax,5)))

    # xticks = [0, 365, 730, 1095, 1460, 1825, 2190]
    # xticklabels = ['2012', '2013', '2014', '2015', '2016', '2017', '2018']
    xticks = [0, 730, 1460, 2190]
    xticklabels = ['2012', '2014', '2016', '2018']

    axes[i, j].set_xticks(xticks)
    axes[i, j].set_xticklabels(xticklabels)
    axes[i, j].spines['top'].set_visible(False)
    axes[i, j].spines['right'].set_visible(False)
    axes[i, j].set_ylabel(f'{units[index]}')

axes[0, 0].legend([line1, line2], ['Simulated data', 'Observed data'], loc='upper left',
                  fontsize=6, frameon=False)

axes[0, 0].text(0.30, 0.05, 'Calibration', color='k', transform=axes[0, 0].transAxes)
axes[0, 0].text(0.70, 0.05, 'Validation', color='k', transform=axes[0, 0].transAxes)

plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.tight_layout()
fig.savefig(f'{all_concentration}\\all_indexes.png', bbox_inches='tight', dpi=500)
fig.savefig(f'{all_concentration}\\all_indexes.pdf', bbox_inches='tight', dpi=500,
            format='pdf')
plt.show()

