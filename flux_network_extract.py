import os.path
import numpy as np
import pandas as pd

from meta import *
from itertools import permutations


def extract_flow_mat(nodes, flux):
    perm = list(permutations(nodes, 2)) + [(x, x) for x in nodes]
    flux_expand = pd.DataFrame(perm, columns=['source', 'target'])
    flux_expand_ = flux_expand.merge(flux, on=['source', 'target'],
                                     how='left')
    flux_expand_ = flux_expand_.fillna(0)
    flux_expand_.columns = ['source', 'target', 'flux_name', 'value']
    flux_mat = flux_expand_.pivot(index='source', columns='target', values='value')
    flux_mat_ = flux_mat.loc[nodes, nodes]
    return flux_mat_


def extract_inout(nodes, flux):
    flux_in = flux.loc[(flux.target.isin(nodes)) & (~flux.source.isin(nodes))]
    flux_out = flux.loc[(flux.source.isin(nodes)) & (~flux.target.isin(nodes))]

    flux_in.columns = ['flux_name', 'source', 'target', 'value']
    flux_in_group = flux_in.loc[:, ['target', 'value']].groupby('target').sum()
    flux_in_group_ = pd.DataFrame(np.zeros(len(nodes)), index=nodes, columns=['value'])
    flux_in_group_.loc[flux_in_group.index] = flux_in_group

    flux_out.columns = ['flux_name', 'source', 'target', 'value']
    flux_out_group = flux_out.loc[:, ['source', 'value']].groupby('source').sum()
    flux_out_group_ = pd.DataFrame(np.zeros(len(nodes)), index=nodes, columns=['value'])
    flux_out_group_.loc[flux_out_group.index] = flux_out_group
    return flux_in_group_, flux_out_group_


def extract_info(mass, flux, scn, save=True, frequency='m'):
    assert frequency in ['d', 'm'], "frequency must be 'd' for daily or " \
                                    "'m' for monthly."
    if frequency == 'm':
        mass = mass.groupby('month').mean()
        flux = flux.groupby('month').mean()

    for t in range(flux.shape[0]):
        flux_Nt = (edges_base.loc[edges_base.type == 'N', ['flux_name', 'source', 'target']]
                   .merge(flux.iloc[t], left_on='flux_name', right_index=True))
        flux_Pt = (edges_base.loc[edges_base.type == 'P', ['flux_name', 'source', 'target']]
                   .merge(flux.iloc[t], left_on='flux_name', right_index=True))

        flux_Nmat = extract_flow_mat(nodes_N, flux_Nt)
        flux_Pmat = extract_flow_mat(nodes_P, flux_Pt)

        input_Nt, output_Nt = extract_inout(nodes_N, flux_Nt)
        input_Pt, output_Pt = extract_inout(nodes_P, flux_Pt)

        mass_Nt = mass.iloc[t].loc[nodes_N]
        mass_Pt = mass.iloc[t].loc[nodes_P]

        if save:
            flux_Nmat.to_csv(f'{network_tables}\\S{scn:02d}\\flux_N_{frequency}{t}.csv')
            flux_Pmat.to_csv(f'{network_tables}\\S{scn:02d}\\flux_P_{frequency}{t}.csv')
            input_Nt.to_csv(f'{network_tables}\\S{scn:02d}\\input_N_{frequency}{t}.csv')
            input_Pt.to_csv(f'{network_tables}\\S{scn:02d}\\input_P_{frequency}{t}.csv')
            output_Nt.to_csv(f'{network_tables}\\S{scn:02d}\\output_N_{frequency}{t}.csv')
            output_Pt.to_csv(f'{network_tables}\\S{scn:02d}\\output_P_{frequency}{t}.csv')
            mass_Nt.to_csv(f'{network_tables}\\S{scn:02d}\\mass_N_{frequency}{t}.csv')
            mass_Pt.to_csv(f'{network_tables}\\S{scn:02d}\\mass_P_{frequency}{t}.csv')


scenarios = range(11)

for scn in scenarios:
    if not os.path.exists(f'{network_tables}\\S{scn:02d}\\'):
        os.mkdir(f'{network_tables}\\S{scn:02d}\\')

nodes_base = pd.read_excel(f'{input_data_path}\\nodes_base.xlsx')
edges_base = pd.read_excel(f'{input_data_path}\\edges_base.xlsx')
nodes_N = nodes_base.loc[nodes_base.type == 'N', 'node'].tolist()
nodes_P = nodes_base.loc[nodes_base.type == 'P', 'node'].tolist()

for scn in scenarios:
    mass_CNP = pd.read_excel(f'{flux_output_path}\\mass_NP_S{scn:02d}.xlsx')
    flux_CNP = pd.read_excel(f'{flux_output_path}\\flux_NP_S{scn:02d}.xlsx')

    extract_info(mass_CNP, flux_CNP, scn, save=True, frequency='m')
    extract_info(mass_CNP, flux_CNP, scn, save=True, frequency='d')
    print(scn)
