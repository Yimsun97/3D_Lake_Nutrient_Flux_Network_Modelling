import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pylab import mpl
import matplotlib

from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from math import sqrt
from math import log
from scipy.spatial import Voronoi, voronoi_plot_2d
from meta import *
from obs_compare import read_result

import geopandas as gpd
import shapely

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


def calc_rmse(y_sim, y_obs):
    y_sim_day = y_sim.copy()
    y_sim_day['TIME'] = y_sim_day['TIME'].astype(int)
    y_obs = y_obs.dropna()
    y_sim_ft = y_sim_day.loc[y_sim_day.TIME.isin(y_obs.DAY.tolist())]
    y_obs_ft = y_obs.loc[y_obs.DAY.isin(y_sim_day.TIME.tolist())]
    return sqrt(mean_squared_error(y_sim_ft.iloc[:, -1].values, y_obs_ft.iloc[:, -1].values))


def create_voronoi(plot=False, saveshp=False):
    dianchi = gpd.read_file('../../概念图/滇池矢量/滇池湖体边界轮廓/Dianchi.shp')
    stat_shp = gpd.read_file('../../概念图/滇池矢量/滇池湖体国控点/湖体国控点84.shp')

    bound = dianchi.geometry[0].buffer(12000).envelope.boundary  # Create a large rectangle surrounding it
    # Create many points along the rectangle boundary. I create one every 100 m.
    boundarypoints = [bound.interpolate(distance=d) for d in range(0, np.ceil(bound.length).astype(int), 100)]
    boundarycoords = np.array([[p.x, p.y] for p in boundarypoints])

    points = np.stack([stat_shp.geometry.x, stat_shp.geometry.y], axis=1)

    all_coords = np.concatenate(
        (boundarycoords, points))  # Create an array of all points on the boundary and inside the polygon

    vor = Voronoi(all_coords)

    lines = [shapely.geometry.LineString(vor.vertices[line]) for line in
             vor.ridge_vertices if -1 not in line]

    polys = shapely.ops.polygonize(lines)
    voronois = gpd.GeoDataFrame(geometry=gpd.GeoSeries(polys))

    poly_gdf = gpd.GeoDataFrame(geometry=[dianchi.geometry[0]])
    points_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(
        x=points[:, 0], y=points[:, 1])
    )

    result = gpd.overlay(df1=voronois, df2=poly_gdf, how="intersection")
    result = result.set_crs(epsg=32648).to_crs(epsg=4326)

    if plot:
        fig, ax = plt.subplots(figsize=(5 / 2.54, 7 / 2.54))
        poly_gdf.boundary.plot(ax=ax, edgecolor="blue", linewidth=2)
        result.plot(ax=ax, color="red", alpha=0.3, edgecolor="black")
        points_gdf.plot(ax=ax, color="maroon")
        plt.show()

    result['AREA'] = result.to_crs(epsg=6933).geometry.area / 1e6
    result['NAME'] = [
        'CGYS', 'LJY', 'WGYS', 'EGYS', 'BYK', 'SDC', 'WHK', 'CHW'
    ]
    if saveshp:
        result.to_file(f'{si_figures}/wq_voronoi.gpkg', driver='GPKG')
    return result


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
for station in stations.keys():
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

# add voronoi
voronois_stats = create_voronoi(plot=False, saveshp=True)
voronois_stats['NameLong'] = voronois_stats['NAME'].map({'CHW': 'Huiwanzhong', 'LJY': 'Luojiaying',
                                                         'EGYS': 'GYSdong', 'CGYS': 'GYSzhong',
                                                         'WGYS': 'GYSxi', 'BYK': 'Baiyukou',
                                                         'SDC': 'Dianchinan', 'WHK': 'Haikouxi'})
voronois_stats['weight'] = voronois_stats['AREA'] / voronois_stats['AREA'].sum()

wq_simu_all = {}
wq_obs_all = {}
for index in indexes:
    wq_simu_list = []
    wq_obs_list = []
    for sta in wq_simu.keys():
        wq_weight = voronois_stats.loc[voronois_stats['NameLong'] == sta, 'weight'].values[0]
        wq_simu_weighted = wq_simu[sta][['TIME', index]].copy()
        wq_simu_weighted[index] = wq_simu_weighted[index] * wq_weight
        wq_simu_list.append(wq_simu_weighted)
        if index != 'ELE':
            wq_obs_weighted = wq_obs[sta][['DAY', index]].copy()
            wq_obs_weighted[index] = wq_obs_weighted[index] * wq_weight
            wq_obs_list.append(wq_obs_weighted)

    wq_simu_sum = np.sum(np.concatenate([x.iloc[:, [-1]].values for x in wq_simu_list], axis=1), axis=1)
    wq_simu_df = pd.DataFrame(wq_simu_sum, columns=[index], index=wq_simu_list[0]['TIME'])
    wq_simu_all[index] = wq_simu_df.reset_index()
    if index != 'ELE':
        wq_obs_sum = np.sum(np.concatenate([x.iloc[:, [-1]].values for x in wq_obs_list], axis=1), axis=1)
        wq_obs_df = pd.DataFrame(wq_obs_sum, columns=[index], index=wq_obs_list[0]['DAY'])
        wq_obs_all[index] = wq_obs_df.reset_index()

fig, axes = plt.subplots(2, 3, figsize=(14 / 2.54, 7 / 2.54), sharex=True)
# plt.tick_params(labelsize=12)

for i, j, index in zip([0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2], indexes):

    axes[i, j].fill_between(wq_simu_all[index]['TIME'],
                            wq_simu_all[index].iloc[:, -1] * (1 - bounds[index]),
                            wq_simu_all[index].iloc[:, -1] * (1 + bounds[index]), alpha=0.5, facecolor=colors[2])
    line1, = axes[i, j].plot(wq_simu_all[index]['TIME'], wq_simu_all[index].iloc[:, -1],
                             c=colors[0], alpha=alpha, label='Simulated data')

    if index == 'ELE':
        line2, = axes[i, j].plot(ELE['DAY'], ELE[index], 'o', c=colors[1],
                                 markersize=2, markeredgewidth=0,
                                 alpha=alpha, label='Observed data')
        # ymax = max((wq_simu_all[index].iloc[:, -1] * (1 + bounds[index])).max(), ELE[index].max())
        ymax = 1890
        axes[i, j].set_ylim(1884, 1890)
    else:
        line2, = axes[i, j].plot(wq_obs_all[index]['DAY'], wq_obs_all[index][index], 'o', c=colors[1],
                                 markersize=2, markeredgewidth=0,
                                 alpha=alpha, label='Observed data')
        ymax = max((wq_simu_all[index].iloc[:, -1] * (1 + bounds[index])).max(), wq_obs_all[index][index].max())
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
fig.savefig(f'{si_figures}\\wq_all_stations.png', bbox_inches='tight', dpi=500)
fig.savefig(f'{si_figures}\\wq_all_stations.pdf', bbox_inches='tight', dpi=500,
            format='pdf')
plt.show()

rmse_dict = {}
for index in indexes:
    if index != 'ELE':
        rmse = calc_rmse(wq_simu_all[index], wq_obs_all[index])
        rmse_dict[index] = rmse
        print(f'{index}: {rmse}')
