import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
from meta import *
import seaborn as sns

sns.set_style("ticks")

font = {'family': 'Arial',
        'size': 7}

matplotlib.rc('font', **font)
plt.rcParams['lines.linewidth'] = 0.5
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

# figure size
figsize = (9 / 2.54, 9 / 2.54)

scn = 0
STAT_period = 'yr_mon'
mass_NP0 = pd.read_excel(f'{flux_output_path}\\mass_NP_S{scn:02d}.xlsx')
flux_NP0 = pd.read_excel(f'{flux_output_path}\\flux_NP_S{scn:02d}.xlsx')

mass_NP0 = mass_NP0.assign(
    TN=mass_NP0.NH4 + mass_NP0.NO3 + mass_NP0.DON + mass_NP0.PON + mass_NP0.ALGAEN,
    TP=mass_NP0.PO4 + mass_NP0.DOP + mass_NP0.POP + mass_NP0.ALGAEP,
    Chla=mass_NP0.ALGAEN,
)
flux_NP0 = flux_NP0.assign(
    TN_PS=flux_NP0.ALGAEN_PS + flux_NP0.DON_PS + flux_NP0.PON_PS + flux_NP0.RNH4_PS + flux_NP0.RNO3_PS,
    TP_PS=flux_NP0.ALGAEP_PS + flux_NP0.DOP_PS + flux_NP0.POP_PS + flux_NP0.PO4_PS,
)
edge_STAT0 = pd.pivot_table(flux_NP0, index=STAT_period, aggfunc='mean').reset_index()
edge_STAT0 = edge_STAT0.assign(
    TN_PS=edge_STAT0.ALGAEN_PS + edge_STAT0.DON_PS + edge_STAT0.PON_PS + edge_STAT0.RNH4_PS + edge_STAT0.RNO3_PS,
    TP_PS=edge_STAT0.ALGAEP_PS + edge_STAT0.DOP_PS + edge_STAT0.POP_PS + edge_STAT0.PO4_PS,
)

month_vis = [5, 6, 7, 8, 9, 10]
month_vis = np.arange(1, 13)
scn_vis = [1, 6, 3, 7]

NP_all = []
stat_all = []
# scn_vis = [6]
for scn in scn_vis:
    mass_NP1 = pd.read_excel(f'{flux_output_path}\\mass_NP_S{scn:02d}.xlsx')
    flux_NP1 = pd.read_excel(f'{flux_output_path}\\flux_NP_S{scn:02d}.xlsx')

    mass_NP1 = mass_NP1.assign(
        TN=mass_NP1.NH4 + mass_NP1.NO3 + mass_NP1.DON + mass_NP1.PON + mass_NP1.ALGAEN,
        TP=mass_NP1.PO4 + mass_NP1.DOP + mass_NP1.POP + mass_NP1.ALGAEP,
        Chla=mass_NP1.ALGAEN,
    )
    mass_NP1.loc[:, 'TN_rate'] = (mass_NP1.TN - mass_NP0.TN) / mass_NP0.TN
    mass_NP1.loc[:, 'TP_rate'] = (mass_NP1.TP - mass_NP0.TP) / mass_NP0.TP
    mass_NP1.loc[:, 'Chla_rate'] = (mass_NP1.Chla - mass_NP0.Chla) / mass_NP0.Chla

    flux_NP1 = flux_NP1.assign(
        TN_PS=flux_NP1.ALGAEN_PS + flux_NP1.DON_PS + flux_NP1.PON_PS + flux_NP1.RNH4_PS + flux_NP1.RNO3_PS,
        TP_PS=flux_NP1.ALGAEP_PS + flux_NP1.DOP_PS + flux_NP1.POP_PS + flux_NP1.PO4_PS,
    )
    flux_NP1.loc[:, 'TN_PS_rate'] = (flux_NP1.TN_PS - flux_NP0.TN_PS) / flux_NP0.TN_PS
    flux_NP1.loc[:, 'TP_PS_rate'] = (flux_NP1.TP_PS - flux_NP0.TP_PS) / flux_NP0.TP_PS

    # aggregate based on STAT_period
    node_STAT1 = pd.pivot_table(mass_NP1, index=STAT_period, aggfunc='mean').reset_index()
    edge_STAT1 = pd.pivot_table(flux_NP1, index=STAT_period, aggfunc='mean').reset_index()

    flux_NP1.loc[:, 'sc'] = f'S{scn:d}'
    edge_STAT1.loc[:, 'sc'] = f'S{scn:d}'

    mass_NP1 = mass_NP1.loc[mass_NP1.month.isin(month_vis)].copy()
    flux_NP1 = flux_NP1.loc[flux_NP1.month.isin(month_vis)].copy()
    NP_all.append(pd.concat([mass_NP1, flux_NP1.iloc[:, 1:]], axis=1))
    stat_all.append(pd.concat([node_STAT1, edge_STAT1.iloc[:, 1:]], axis=1))

NP_df = pd.concat(NP_all, axis=0, ignore_index=True)
NP_df_ = NP_df.loc[:, ['TN_PS_rate', 'TP_PS_rate', 'TN_rate', 'TP_rate', 'Chla_rate', 'sc']].copy()
NP_melted = pd.melt(NP_df_, id_vars='sc', value_vars=['TN_PS_rate', 'TP_PS_rate', 'TN_rate', 'TP_rate',
                                                      'Chla_rate'])

NP_stat_df = pd.concat(stat_all, axis=0, ignore_index=True)
NP_stat_df_ = NP_stat_df.loc[:, ['TN_PS_rate', 'TP_PS_rate', 'TN_rate', 'TP_rate', 'Chla_rate', 'sc']].copy()
NP_stat_melted = pd.melt(NP_stat_df_, id_vars='sc', value_vars=['TN_PS_rate', 'TP_PS_rate', 'TN_rate', 'TP_rate',
                                                                'Chla_rate'])

# network indicators
network_inds = pd.read_csv('ENA/NetworkIndicators.csv', header=0)
network_rtst = network_inds.loc[network_inds['index'] == 'rTST'].copy()
s0_rtst_N = network_rtst.loc[(network_rtst.sc == 'S00') & (network_rtst.varname == 'N'),
                             'value'].values
s0_rtst_P = network_rtst.loc[(network_rtst.sc == 'S00') & (network_rtst.varname == 'P'),
                             'value'].values
ratio_df_list = []
for scn in scn_vis:
    sc_rtst_N = network_rtst.loc[(network_rtst.sc == f'S0{scn}') & (network_rtst.varname == 'N'),
                                 'value'].values
    sc_rtst_P = network_rtst.loc[(network_rtst.sc == f'S0{scn}') & (network_rtst.varname == 'P'),
                                 'value'].values
    ratio_N = (sc_rtst_N - s0_rtst_N) / s0_rtst_N
    ratio_P = (sc_rtst_P - s0_rtst_P) / s0_rtst_P

    ratio_N_df = pd.DataFrame(ratio_N, columns=['value'])
    ratio_P_df = pd.DataFrame(ratio_P, columns=['value'])

    ratio_N_df.insert(0, 'variable', 'rTST of N cycle')
    ratio_P_df.insert(0, 'variable', 'rTST of P cycle')

    ratio_N_df.insert(0, 'sc', f'S{scn}')
    ratio_P_df.insert(0, 'sc', f'S{scn}')

    ratio_df = pd.concat([ratio_N_df, ratio_P_df], axis=0, ignore_index=True)
    ratio_df_list.append(ratio_df)

ratio_df_all = pd.concat(ratio_df_list, axis=0, ignore_index=True)

# combine NP_melted and ratio_df_all
# NP_melted = pd.concat([NP_melted, ratio_df_all], axis=0, ignore_index=True)

# convert to percentage %
NP_melted.loc[:, 'value'] = NP_melted.value * 100
NP_stat_melted.loc[:, 'value'] = NP_stat_melted.value * 100


# convert legend names
varname_dict = {
    'TN_PS_rate': 'N inflow', 'TP_PS_rate': 'P inflow',
    'TN_rate': 'TN', 'TP_rate': 'TP',
    'Chla_rate': 'Chla',
}

for key, value in varname_dict.items():
    NP_melted.loc[NP_melted.variable == key, 'variable'] = value
    NP_stat_melted.loc[NP_stat_melted.variable == key, 'variable'] = value

# change colors
color_dict = {
    'N inflow': '#2171b5', 'P inflow': '#e6754d',
    'TN': '#72a7ce', 'TP': '#ecb98e',
    'Chla': '#756bb1',
}
palette = color_dict
# palette = "dark"

# # Draw a nested barplot
# g = sns.catplot(
#     data=NP_melted, kind="bar",
#     x="variable", y="value", hue="sc", capsize=.05,
#     errorbar="sd", palette="dark", alpha=.6, height=6
# )
# # g.despine(left=True)
# g.set_axis_labels("Scenario", "Change relative to baseline (%)")
# g.legend.set_title("")
# plt.show()

fig, ax = plt.subplots(figsize=figsize, dpi=500)
g = sns.barplot(
    ax=ax,
    data=NP_melted,
    x="sc", y="value", hue="variable", capsize=.03,
    errorbar="sd", palette=palette, alpha=.8,
)
sns.move_legend(g, "center right", bbox_to_anchor=(0.8, 0.8))
ax.set_xlabel("Scenario")
ax.set_ylabel("Change relative to baseline (%)")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.legend(frameon=False)
plt.setp(ax.patches, linewidth=0)
plt.savefig(f'{all_concentration}\\nonlinearity_1367.pdf', dpi=500, format='pdf')
plt.show()
