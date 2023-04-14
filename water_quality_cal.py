import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
import pylab as mpl
from meta import *

mpl.rcParams['font.sans-serif'] = ['SimSun']  # 字体 SimSun，FangSong，SimHei，Times New Roman
mpl.rcParams['axes.unicode_minus'] = False  # 解决保存图像时'-'显示为方块的问题
plt.rcParams['font.size'] = 8  # 单位为磅
plt.rcParams['lines.linewidth'] = 1  # 线宽
plt.rcParams['lines.linestyle'] = '-'  # 线型   - -- -. :
plt.rcParams['lines.markersize'] = 3  # 标记点大小
plt.rcParams['lines.marker'] = 'None'  # 标记类型 None . , o v ^ < > s + x D d
plt.rcParams['mathtext.fontset'] = 'stix'  # 公式中的字体，stix与新罗马字体很接近
# mpl.rcParams.update(mpl.rcParamsDefault) #全部恢复默认值

scenarios = range(11)  # 'S3','S3a'

station = stations[0]
ylabels = {'NH4': r'$\mathregular{NH_3-N}$', 'TN': 'TN', 'TP': 'TP', 'DO': 'DO', 'CHLA': r'Chla'}
results = {'TN': pd.DataFrame(), 'TP': pd.DataFrame(), 'NH4': pd.DataFrame(), 'DO': pd.DataFrame(),
           'CHLA': pd.DataFrame(), 'ELE': pd.DataFrame()}
data_range = {'TN': pd.DataFrame(), 'TP': pd.DataFrame(), 'NH4': pd.DataFrame(), 'DO': pd.DataFrame(),
              'CHLA': pd.DataFrame(), 'ELE': pd.DataFrame()}
wqs = ['TN', 'CHLA', 'TP', 'NH4', 'ELE']
new_names = {'ALGAEN_MASS': 'ALGAEN', 'PON_MASS': 'PON', 'DON_MASS': 'DON',
             'NH4_MASS': 'NH4', 'NO3_MASS': 'NO3',
             'ALGAEP_MASS': 'ALGAEP', 'PO4_MASS': 'PO4', 'TPOP_MASS': 'POP',
             'TDOP_MASS': 'DOP', 'LDP_FROM_ALGAE': 'DOP_FROM_ALGAE',
             'P_ALGAE_UPTK': 'PO4_ALGAE_UPTK', 'DOP_FROM_LOP': 'DOP_FROM_POP'}

for scn in scenarios:
    wqwcts = pd.read_table(f'{scn_path}\\S{scn:02d}\\WQWCTS.OUT', delim_whitespace=True)
    wqwcts = wqwcts[wqwcts['K'] == wqwcts['K'].max()]
    wqwcts = wqwcts.loc[(wqwcts['TIME'] % 1 > 0.45) & (wqwcts['TIME'] % 1 < 0.55) &
                        (wqwcts['I'] == station[0]) & (wqwcts['J'] == station[1])]
    wqwcts.index = range(len(wqwcts))

    N = pd.read_table(f'{scn_path}\\S{scn:02d}\\WQ_NBALANCE.OUT', sep='\s+')
    P = pd.read_table(f'{scn_path}\\S{scn:02d}\\WQ_PBALANCE1.OUT', sep='\s+')
    flux_inp = pd.concat([N, P.iloc[:, 1:]], axis=1)
    flux_inp = flux_inp.rename(columns=new_names)
    mass = flux_inp[['JDAY', 'ALGAEN', 'ALGAEP', 'PO4', 'DOP', 'POP', 'NH4', 'NO3', 'DON', 'PON']]
    mass = mass[mass['JDAY'] % 1 > 0.45]
    mass.index = range(len(mass))
    mass['ELE'] = wqwcts['ELE']
    # 根据dxdy.inp计算
    mass['VL'] = 1487500379 + (mass['ELE'] - 1887.4531) * 295547309.7
    data = mass.copy()
    data['CHLA'] = 1000000 * mass['ALGAEN'] / 0.16 / 0.025 / 1000 / mass['VL']
    data['NH4'] = 1000000 * mass['NH4'] / mass['VL']
    data['TN'] = 1000000 * mass[['ALGAEN', 'NH4', 'NO3', 'PON', 'DON']].sum(axis=1) / mass['VL']
    data['TP'] = 1000000 * mass[['ALGAEP', 'PO4', 'DOP', 'POP']].sum(axis=1) / mass['VL']

    for wq in wqs:
        results[wq] = data[wq]

    for wq in wqs:
        data_range[wq]['TIME'] = data['JDAY']
        data_range[wq]['%s' % scn] = results[wq]
    print(scn)

for wq in wqs:
    data_range[wq]['date'] = data_range[wq]['TIME'].apply(lambda x: pd.Timedelta(days=x - 1) +
                                                                    pd.to_datetime(jday_start))
    data_range[wq]['month'] = data_range[wq]['date'].apply(lambda x: x.month)
    monthly = pd.pivot_table(data_range[wq], index='month', aggfunc='mean').reset_index()
    data_range[wq].to_excel(f'{wq_output_path}\\wq_max+min_{wq}.xlsx', index=False)
    monthly.to_excel(f'{wq_output_path}\\monthly_{wq}.xlsx', index=False)
