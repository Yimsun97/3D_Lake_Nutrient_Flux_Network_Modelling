import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
import pylab as mpl
import matplotlib.dates as mdates
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


scenarios = [x for x in range(11)]
scenarios_fig4 = [0, 1, 6]
scenarios_fig5 = [0, 3, 7]
scenarios_fig_si1 = [0, 2, 4, 5]
scenarios_fig_si2 = [0, 8, 9, 10]
labels_fig4 = [f'S{x}' for x in scenarios_fig4]
labels_fig5 = [f'S{x}' for x in scenarios_fig5]
labels_fig_si1 = [f'S{x}' for x in scenarios_fig_si1]
labels_fig_si2 = [f'S{x}' for x in scenarios_fig_si2]

clr = list(plt.get_cmap('tab10').colors)
# clr = ['#ad1d1f', '#eb252d', '#f8ae1f', '#0a467b', '#00b0db', '#705491']
clr.insert(0, (0.0, 0.0, 0.0))

cmap = plt.get_cmap('tab10', 10)
hexs = []
for i in range(cmap.N):
    rgba = cmap(i)
    # rgb2hex accepts rgb or rgba
    hexs.append(matplotlib.colors.rgb2hex(rgba))

ylabels = {'TEM': 'Temperature (°C)', 'COD': '$COD\ (mg\ L^{-1})$',
           'NH4': r'$NH_3-N\ (mg\ L^{-1})$', 'ELE': '水位 (m)',
           'TN': '$TN\ (mg\ L^{-1})$', 'TP': '$TP\ (mg\ L^{-1})$',
           'DO': '$DO\ (mg\ L^{-1})$', 'CHLA': r'$Chla\ (mg\ L^{-1})$'}
data_range = {}
wqs = ['TN', 'TP', 'CHLA']
for wq in wqs:
    data = pd.read_excel(f'{wq_output_path}\\wq_max+min_{wq}.xlsx')
    data_range[wq] = data

# %% Fig. 4a-4c
fig, axes = plt.subplots(1, len(wqs), figsize=(19 / 2.54, 6 / 2.54),
                         dpi=500)

for j, wq in enumerate(wqs):

    sc_ids = scenarios_fig4
    labels = labels_fig4

    colors = [clr[x] for x in sc_ids]
    for i in range(len(sc_ids)):
        axes[j].plot(data_range[wq]['date'], data_range[wq]['%s' % sc_ids[i]],
                     label=labels[i], c=colors[i], alpha=1.0)

    axes[j].yaxis.set_major_formatter(FormatStrFormatter('%3.2f'))
    axes[j].set_ylabel(ylabels[wq])
    axes[j].set(xlabel=None)
    axes[j].tick_params(axis='x', labelrotation=30)
    axes[j].spines['top'].set_visible(False)
    axes[j].spines['right'].set_visible(False)


axes[-1].legend(loc=2, frameon=False)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))

text_x = 0.9
text_y = 0.9
axes[0].text(text_x, text_y, '(a)', transform=axes[0].transAxes)
axes[1].text(text_x, text_y, '(b)', transform=axes[1].transAxes)
axes[2].text(text_x, text_y, '(c)', transform=axes[2].transAxes)

fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
fig.tight_layout()
fig.savefig(f'{wq_line_plot_path}\\WQ_sc_fig4.jpg', dpi=500)
fig.savefig(f'{wq_line_plot_path}\\WQ_sc_fig4.pdf', dpi=500, format='pdf')

plt.show()

# %% Fig. 5a-5c
fig, axes = plt.subplots(1, len(wqs), figsize=(19 / 2.54, 6 / 2.54),
                         dpi=500)

for j, wq in enumerate(wqs):

    sc_ids = scenarios_fig5
    labels = labels_fig5

    colors = [clr[x] for x in sc_ids]
    for i in range(len(sc_ids)):
        axes[j].plot(data_range[wq]['date'], data_range[wq]['%s' % sc_ids[i]],
                     label=labels[i], c=colors[i], alpha=1.0)

    axes[j].yaxis.set_major_formatter(FormatStrFormatter('%3.2f'))
    axes[j].set_ylabel(ylabels[wq])
    axes[j].set(xlabel=None)
    axes[j].tick_params(axis='x', labelrotation=30)
    axes[j].spines['top'].set_visible(False)
    axes[j].spines['right'].set_visible(False)


axes[-1].legend(loc=2, frameon=False)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))

text_x = 0.9
text_y = 0.9
axes[0].text(text_x, text_y, '(a)', transform=axes[0].transAxes)
axes[1].text(text_x, text_y, '(b)', transform=axes[1].transAxes)
axes[2].text(text_x, text_y, '(c)', transform=axes[2].transAxes)

fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
fig.tight_layout()
fig.savefig(f'{wq_line_plot_path}\\WQ_sc_fig5.jpg', dpi=500)
fig.savefig(f'{wq_line_plot_path}\\WQ_sc_fig5.pdf', dpi=500, format='pdf')

plt.show()


# %% Fig. si 1
fig, axes = plt.subplots(1, len(wqs), figsize=(19 / 2.54, 6 / 2.54),
                         dpi=500)

for j, wq in enumerate(wqs):

    sc_ids = scenarios_fig_si1
    labels = labels_fig_si1

    colors = [clr[x] for x in sc_ids]
    for i in range(len(sc_ids)):
        axes[j].plot(data_range[wq]['date'], data_range[wq]['%s' % sc_ids[i]],
                     label=labels[i], c=colors[i], alpha=1.0)

    axes[j].yaxis.set_major_formatter(FormatStrFormatter('%3.2f'))
    axes[j].set_ylabel(ylabels[wq])
    axes[j].set(xlabel=None)
    axes[j].tick_params(axis='x', labelrotation=30)
    axes[j].spines['top'].set_visible(False)
    axes[j].spines['right'].set_visible(False)


axes[-1].legend(loc=2, frameon=False)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))

text_x = 0.9
text_y = 0.9
axes[0].text(text_x, text_y, '(a)', transform=axes[0].transAxes)
axes[1].text(text_x, text_y, '(b)', transform=axes[1].transAxes)
axes[2].text(text_x, text_y, '(c)', transform=axes[2].transAxes)

fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
fig.tight_layout()
fig.savefig(f'{wq_line_plot_path}\\WQ_sc_fig_si1.jpg', dpi=500)
fig.savefig(f'{wq_line_plot_path}\\WQ_sc_fig_si1.pdf', dpi=500, format='pdf')

plt.show()

# %% Fig. si 2
fig, axes = plt.subplots(1, len(wqs), figsize=(19 / 2.54, 6 / 2.54),
                         dpi=500)

for j, wq in enumerate(wqs):

    sc_ids = scenarios_fig_si2
    labels = labels_fig_si2

    colors = [clr[x] for x in sc_ids]
    for i in range(len(sc_ids)):
        axes[j].plot(data_range[wq]['date'], data_range[wq]['%s' % sc_ids[i]],
                     label=labels[i], c=colors[i], alpha=1.0)

    axes[j].yaxis.set_major_formatter(FormatStrFormatter('%3.2f'))
    axes[j].set_ylabel(ylabels[wq])
    axes[j].set(xlabel=None)
    axes[j].tick_params(axis='x', labelrotation=30)
    axes[j].spines['top'].set_visible(False)
    axes[j].spines['right'].set_visible(False)


axes[-1].legend(loc=2, frameon=False)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))

text_x = 0.9
text_y = 0.9
axes[0].text(text_x, text_y, '(a)', transform=axes[0].transAxes)
axes[1].text(text_x, text_y, '(b)', transform=axes[1].transAxes)
axes[2].text(text_x, text_y, '(c)', transform=axes[2].transAxes)

fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
fig.tight_layout()
fig.savefig(f'{wq_line_plot_path}\\WQ_sc_fig_si2.jpg', dpi=500)
fig.savefig(f'{wq_line_plot_path}\\WQ_sc_fig_si2.pdf', dpi=500, format='pdf')

plt.show()
