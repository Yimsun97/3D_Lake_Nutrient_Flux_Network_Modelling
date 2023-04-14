import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
from meta import *

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

STAT_period = 'yr_mon'

# scale factor for N and P for visualization
scale_N_MASS = 1.2
scale_P_MASS = 13.2
scale_N_FLUX = 0.3
scale_P_FLUX = 5.0

# node variables
node_variables = ['ALGAEN', 'ALGAEP', 'PO4', 'DOP', 'POP', 'NH4', 'NO3', 'DON', 'PON',
                  'SED_P', 'SED_N']

# label language and font family
label_id = 'label0'
unit_text = 'Unit: Tons'

# figure size
figsize = (19 / 2 / 2.54, 19 / 2 / 2.54)
# font size
fontsize = 7
font_size = 7
# figure show on or off
fig_show = True
# alpha
alpha_node = 0.8
alpha_edge = 0.6

scn = 0
STAT_period = 'yr_mon'
mass_NP0 = pd.read_excel(f'{flux_output_path}\\mass_NP_S{scn:02d}.xlsx')
flux_NP0 = pd.read_excel(f'{flux_output_path}\\flux_NP_S{scn:02d}.xlsx')
# flux_NP0[flux_NP0.columns[1:-7]] = flux_NP0[flux_NP0.columns[1:-7]].abs()

nodes_base = pd.read_excel(f'{input_data_path}\\nodes_base.xlsx')
nodes_base.index = nodes_base['node']
nodes_base['color'] = '#aaaaaa'
edges_base = pd.read_excel(f'{input_data_path}\\edges_base.xlsx')
edges_base.index = edges_base['flux_name']
node_N_names = [i for i in nodes_base.index if nodes_base['type'][i] == 'N']
node_P_names = [i for i in nodes_base.index if nodes_base['type'][i] == 'P']
edge_N_names = [i for i in edges_base.index if edges_base['type'][i] == 'N']
edge_P_names = [i for i in edges_base.index if edges_base['type'][i] == 'P']

scale_N_MASS = 10
scale_P_MASS = scale_N_MASS * 8
scale_N_FLUX = 1
scale_P_FLUX = scale_N_FLUX * 8

scn_vis = np.arange(11)
# scn_vis = [6]
for scn in scn_vis:
    mass_NP1 = pd.read_excel(f'{flux_output_path}\\mass_NP_S{scn:02d}.xlsx')
    flux_NP1 = pd.read_excel(f'{flux_output_path}\\flux_NP_S{scn:02d}.xlsx')
    # flux_NP1[flux_NP1.columns[1:-7]] = flux_NP1[flux_NP1.columns[1:-7]].abs()
    mass_NP = mass_NP1.iloc[:, :-7] - mass_NP0.iloc[:, :-7]
    mass_NP[mass_NP1.iloc[:, -7:].columns] = mass_NP1.iloc[:, -7:]
    flux_NP = flux_NP1.iloc[:, :-7] - flux_NP0.iloc[:, :-7]
    flux_NP[flux_NP1.iloc[:, -7:].columns] = flux_NP1.iloc[:, -7:]

    # aggregate based on STAT_period
    node_STAT = pd.pivot_table(mass_NP, index=STAT_period,
                               values=node_variables,
                               aggfunc='mean').reset_index()
    edge_STAT = pd.pivot_table(flux_NP, index=STAT_period, aggfunc='mean').reset_index()
    edge_STAT1 = pd.pivot_table(flux_NP1, index=STAT_period, aggfunc='mean').reset_index()

    # scale the nodes and edges for visualization
    node_times = node_STAT.copy()
    edge_times = edge_STAT.copy()
    node_times.loc[:, node_N_names] = node_times[node_N_names] * scale_N_MASS
    node_times.loc[:, node_P_names] = node_times[node_P_names] * scale_P_MASS
    edge_times.loc[:, edge_N_names] = edge_times[edge_N_names] * scale_N_FLUX
    edge_times.loc[:, edge_P_names] = edge_times[edge_P_names] * scale_P_FLUX

    # position of stock nodes and labels
    pos = dict(zip(nodes_base['node'], nodes_base[['posX', 'posY']].values))

    # position of boundary nodes and labels
    nodes_base['posY1'] = nodes_base['posY'] + 0.05
    pos_boundary = dict(zip(nodes_base['node'], nodes_base[['posX', 'posY1']].values))

    t = [x - 1 for x in wet_seasons]
    # update node sizes and labels
    # Note: node sizes are scaled proportional to true values
    nodes_season_mean = node_times.loc[t, :].iloc[:, 1:].mean(axis=0)
    nodes_STAT_season_mean = node_STAT.loc[t, :].iloc[:, 1:].mean(axis=0)
    nodes_base['size'].update(nodes_season_mean)
    nodes_base['size'] = nodes_base['size'].clip(*nodes_clip)

    for i in node_times.columns[1:]:
        if nodes_season_mean.loc[i] > 0:
            nodes_base.loc[i, 'color'] = '#c12e34'
        elif nodes_season_mean.loc[i] < 0:
            nodes_base.loc[i, 'color'] = '#0098d9'
        else:
            nodes_base.loc[i, 'color'] = 'grey'

    # update edge widths and labels
    # Note: edge widths are scaled proportional to true values
    edges_season_mean = edge_times.loc[t, :].iloc[:, 1:-1].abs().mean(axis=0)
    edges_base.loc[:, 'size'].update(edges_season_mean)
    edges_base.loc[:, 'size'] = edges_base['size'].clip(*edges_clip)
    edges_base.loc[:, 'label'] = edges_season_mean.apply(lambda x: f'{x:.1f}')

    # change arrow directions if necessary
    edges_base_ = edges_base.copy()
    changed_fluxes = ['RNH4_BEN1', 'RNO3_BEN1', 'PO4_BEN1']

    for flux_ in changed_fluxes:
        if edge_STAT1.loc[t, flux_].mean() > 0:
            src = edges_base_.loc[flux_, 'source']
            tar = edges_base_.loc[flux_, 'target']
            edges_base_.loc[flux_, 'source'] = tar
            edges_base_.loc[flux_, 'target'] = src

    for i in edges_base_.index:
        if edge_times.loc[t, i].mean() == 0:
            edges_base_.loc[i, 'color'] = 'grey'
        else:
            if edge_times.loc[t, i].mean() > 0:
                if edge_STAT1.loc[t, i].mean() >= 0:
                    edges_base_.loc[i, 'color'] = '#c12e34'
                else:
                    edges_base_.loc[i, 'color'] = '#0098d9'
            else:
                if edge_STAT1.loc[t, i].mean() >= 0:
                    edges_base_.loc[i, 'color'] = '#0098d9'
                else:
                    edges_base_.loc[i, 'color'] = '#c12e34'

    G = nx.from_pandas_edgelist(edges_base_, source='source', target='target', edge_attr=['size'],
                                create_using=nx.DiGraph())
    edges_base_1 = edges_base_.drop_duplicates(['source', 'target'])
    edges_base_2 = edges_base_.loc[edges_base_.duplicated(['source', 'target'])]
    G1 = nx.from_pandas_edgelist(edges_base_1, source='source', target='target', edge_attr=['size'],
                                 create_using=nx.DiGraph())
    G2 = nx.from_pandas_edgelist(edges_base_2, source='source', target='target', edge_attr=['size'],
                                 create_using=nx.DiGraph())
    # draw stock nodes
    WQnode = nodes_base.dropna(subset=['type']).copy()
    WQnode.loc[:, 'label_size'] = nodes_STAT_season_mean.loc[WQnode.index]
    WQnode.loc[:, 'label'] = WQnode.apply(lambda x: f"{x[label_id]}\n{x['label_size']:.1f}", axis=1)

    fig, axes = plt.subplots(1, 1, figsize=figsize, dpi=500)

    nx.draw_networkx_nodes(G, pos=pos, ax=axes,
                           nodelist=WQnode['node'].values,
                           node_shape='o',
                           node_color=WQnode['color'].values,
                           node_size=WQnode['size'].astype(float).values,
                           alpha=alpha_node,
                           label='Mass in water column',
                           )
    nx.draw_networkx_labels(G, pos=pos, ax=axes,
                            labels=dict(zip(WQnode['node'], WQnode['label'])),
                            font_size=font_size,
                            )

    # draw boundary nodes
    BOUNDARYnode = nodes_base.loc[nodes_base.type.isna()].copy()
    nx.draw_networkx_nodes(G, pos=pos, ax=axes,
                           nodelist=BOUNDARYnode['node'].values,
                           node_shape='s',
                           node_color=BOUNDARYnode['color'].values,
                           node_size=BOUNDARYnode['size'].astype(float).values,
                           alpha=alpha_node,
                           label='Boundary',
                           )
    nx.draw_networkx_labels(G, pos=pos_boundary, ax=axes,
                            labels=dict(zip(BOUNDARYnode['node'], BOUNDARYnode[label_id])),
                            font_size=fontsize,
                            )


    nx.draw_networkx_edges(G1, pos=pos, ax=axes,
                           edgelist=edges_base_1.loc[:, ['source', 'target']].values.tolist(),
                           width=edges_base_1['size'].values,
                           edge_color=edges_base_1['color'].values,
                           node_size=nodes_base['size'].astype(float).values,
                           nodelist=list(nodes_base['node'].values),
                           alpha=alpha_edge,
                           arrows=True,
                           arrowsize=10,
                           arrowstyle='fancy',
                           connectionstyle='arc3, rad=-0.2',
                           label='Flux',
                           )
    nx.draw_networkx_edges(G2, pos=pos, ax=axes,
                           edgelist=edges_base_2.loc[:, ['source', 'target']].values.tolist(),
                           width=edges_base_2['size'].values,
                           edge_color=edges_base_2['color'].values,
                           node_size=nodes_base['size'].astype(float).values,
                           nodelist=list(nodes_base['node'].values),
                           alpha=alpha_edge,
                           arrows=True,
                           arrowsize=10,
                           arrowstyle='fancy',
                           connectionstyle='arc3, rad=0.2',
                           label='Flux',
                           )

    axes.spines['top'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['left'].set_visible(False)

    fig.subplots_adjust(left=0.0, bottom=0.0, right=1, top=1, wspace=0, hspace=0)
    fig.tight_layout()
    fig.savefig(f'{flux_network_diff_path}\\S0S{scn}_wet.jpg', dpi=500)

    if scn in (1, 6, 3, 7):
        fig.savefig(f'{flux_network_diff_path}\\S0S{scn}_wet.pdf', dpi=500, format='pdf')

    if fig_show:
        plt.show()
