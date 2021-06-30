"""Download format data and make a plot"""
import os
import datetime
import pandas as pd
import numpy as np
import seaborn as sns
from multiprocessing import Pool
from subprocess import Popen, PIPE
from io import StringIO
from functools import reduce
from matplotlib import pyplot as plt
from sklearn.preprocessing import label_binarize
from sklearn import svm, datasets
from sklearn.metrics import confusion_matrix
from time import sleep
from sklearn.preprocessing import MinMaxScaler
from decimal import Decimal
from matplotlib.lines import Line2D

pipeline_dir = "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/vogdb_by_vogdb/"
time_paths = [
    os.path.join(
        pipeline_dir,
        "hmmer/vogdb_by_vogdb_hmmer_results_2021_06_18_09/run_times/"),
    # os.path.join(
    #     pipeline_dir,
    #     "mmseqs/vogdb_by_vogdb_mmseqs_t2_results_2021_06_23_11/run_times"),
    os.path.join(
        pipeline_dir,
        "mmseqs/vogdb_by_vogdb_mmseqs_t32_results_2021_06_17_18/run_times/")
]
os.system('/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/vogdb_by_vogdb/hmmer/vogdb_by_vogdb_hmmer_results_2021_06_18_09/run_times/time_hmmer')
reduce(lambda x, y: x + y, )
# FIX NOT FINAL [-1] FOR BAD HMMER STUFF
p  = [os.path.join(t_path, t_file)
 for t_path in time_paths
 for t_file in os.listdir(t_path)]

pd.read_csv('/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/vogdb_by_vogdb/hmmer/vogdb_by_vogdb_hmmer_results_2021_06_18_09/run_times/time_hmmer',
                          sep=": ",
                          header=None,
                          index_col=0).T
for i in p:
    print(i)
    pd.read_csv(i,
                              sep=": ",
                              header=None,
                              index_col=0).T

data = pd.concat([pd.read_csv(os.path.join(t_path, t_file),
                              sep=": ",
                              header=None,
                              index_col=0).T
                  for t_path in time_paths
                  for t_file in os.listdir(t_path)]).reset_index(drop=True)
data['Elapsed (wall clock) time (Seconds)'] = \
    data['Elapsed (wall clock) time (h:mm:ss or m:ss)'].apply(
        lambda x: sum([float(i) * (60**e)
                       for e, i in enumerate(reversed(x.split(':')))])
    )
data['Annotation Tool'] = [
    "MMseqs2 with 32 Threads"
    if i.find('--threads 32') > 0 else "MMseqs2 with 2 Threads"
    if i.find('mmseqs search') > 0 else "HMMER3"
    for i in data['Command being timed'].values]

data['Annotation Algorithm'] = [
    "MMseqs2 with 32 Threads, %s Sensitivity"
    %(lambda x: x[x.index('-s') + 1])(i.split(' '))
    if i.find('--threads 32') > 0 else
    "MMseqs2 with 2 Threads, %s Sensitivity"
    %(lambda x: x[x.index('-s') + 1])(i.split(' '))
    if i.find('mmseqs search') > 0 else "HMMER3"
    for i in data['Command being timed'].values]

data.sort_values('Annotation Algorithm', inplace=True)

data['Percent of CPU this job got'] = \
    data['Percent of CPU this job got'].str.replace('%', '').\
    astype(float)
data['Maximum resident set size (kbytes)'] = \
    data['Maximum resident set size (kbytes)'].astype(float)

mmseqs_legend = np.unique([
    "%s Sensitivity"
    %(lambda x: x[x.index('-s') + 1])(i.split(' '))
    for i in data['Command being timed'].values
    if i.find('mmseqs search') > 0])
sum


lumins = int(sum(data['Command being timed'].str.find('mmseqs search') > 0)\
             / 2)

###############
#  Plot time  #
###############

def make_plot_row(axs, y_val, y_log_scale=False):

    sns.set_palette(sns.color_palette("Reds",3)[2:])
    sns.barplot(y=y_val,
                x='Annotation Tool',
                hue='Annotation Algorithm',
                ax=axs[0],
                data=data[data['Annotation Tool'] == 'HMMER3' ],
                saturation=0.5
                )
    sns.set_palette(sns.color_palette("Greens", lumins + 5)[5:])
    sns.barplot(y=y_val,
                x='Annotation Tool',
                hue='Annotation Algorithm',
                ax=axs[1],
                data=data[data['Annotation Tool'] == 'MMseqs2 with 2 Threads' ],
                saturation=0.4
                )
    sns.set_palette(sns.color_palette("Blues", lumins + 5)[5:])
    sns.barplot(y=y_val,
                x='Annotation Tool',
                hue='Annotation Algorithm',
                ax=axs[2],
                data=data[data['Annotation Tool'] == 'MMseqs2 with 32 Threads' ],
                saturation=0.4
                )

    for i in range(3):
        axs[i].set_xlabel('')
        if y_log_scale:
            axs[i].set_yscale('log')
            axs[i].set_ylim([1, data[y_val].max()])
        else:
            axs[i].set_ylim([0, data[y_val].max()])
        # axs[i].patch.set_facecolor(DARKER_COLOR)

    axs[1].yaxis.set_ticklabels([])
    axs[1].set_ylabel('')
    axs[2].yaxis.set_ticklabels([])
    axs[2].set_ylabel('')
    legend_elements = [
        Line2D([0], [0], color=j, lw=2, label=i)
        for i, j in zip(mmseqs_legend,
                        sns.color_palette("Greys", lumins + 5)[5:])
    ] + [Line2D([0], [0], color=sns.color_palette("Reds",3, desat=0.5)[2],
                lw=2, label="HMMER3")]
    axs[1].legend(handles=legend_elements, loc='upper center', ncol=6,
                  bbox_to_anchor=(0.9, -0.1))
    axs[0].legend([],[], frameon=False)
    axs[2].legend([],[], frameon=False)
    return axs


def make_plot_just_one(ax, y_val):
    """
    :param ax:
    :param y_val:
    :returns:
    """
    sns.set_palette(sns.color_palette("Blues", lumins + 5)[5:])
    sns.barplot(y=y_val,
                x='Annotation Tool',
                hue='Annotation Algorithm',
                ax=ax,
                data=data[data['Annotation Tool'] == 'MMseqs2 with 32 Threads' ],
                saturation=0.4
                )

    ax.set_xlabel('')
    legend_elements = [
        Line2D([0], [0], color=j, lw=2, label=i)
        for i, j in zip(mmseqs_legend,
                        sns.color_palette("Blues", lumins + 5)[5:])
    ]
    ax.legend(handles=legend_elements, loc='upper center', ncol=3,
                  bbox_to_anchor=(0.45, -0.1))
    # ax.patch.set_facecolor(DARKER_COLOR)
    return ax


###############
#  Plot Clock #
###############

sns.set_style("whitegrid")
    #sns.despine(bottom = True, left = True)

plt.clf()
fig, axs = plt.subplots(figsize=(12, 6), ncols=3,
                        gridspec_kw={'width_ratios': [2, lumins, lumins]})
axs = make_plot_row(axs, 'Elapsed (wall clock) time (Seconds)',
                    y_log_scale=True)
fig.patch.set_facecolor(LIGHT_COLOR)
fig.subplots_adjust(wspace=0, bottom=0.25)
fig.suptitle("Computation Speed of MMseqs2 and HMMER3, with E-Value = 1e-15",
             fontsize=16)
fig.savefig('plots/speed_results_clock.png')

###############
#  Plot CPU % #
###############


plt.clf()
fig, axs = plt.subplots(figsize=(12, 6), ncols=3,
                        gridspec_kw={'width_ratios': [2, lumins, lumins]})
axs = make_plot_row(axs, 'Percent of CPU this job got')
fig.patch.set_facecolor(LIGHT_COLOR)
fig.subplots_adjust(wspace=0, bottom=0.25)
fig.suptitle("CPU performance of MMseqs2 and hmmer3, with e-value = 1e-15",
             fontsize=16)
fig.savefig('plots/speed_results_cpu.png')

#################
#  Plot Memory  #
#################

plt.clf()
fig, axs = plt.subplots(figsize=(12, 6), ncols=3,
                        gridspec_kw={'width_ratios': [2, lumins, lumins]})
axs = make_plot_row(axs, 'Maximum resident set size (kbytes)')
fig.patch.set_facecolor(LIGHT_COLOR)
fig.subplots_adjust(wspace=0, bottom=0.25)
fig.suptitle("Memory Uses of MMseqs2 and HMMER3, with E-Value = 1e-15",
             fontsize=16)
fig.savefig('plots/speed_results_memory.png')

data.set_index('Annotation Algorithm', inplace=True)

seconds32 = data.loc['HMMER3', 'Elapsed (wall clock) time (Seconds)'] -\
    data.loc['MMseqs2 with 32 Threads, 7.5 Sensitivity',
             'Elapsed (wall clock) time (Seconds)']
seconds2 = data.loc['HMMER3', 'Elapsed (wall clock) time (Seconds)'] -\
    data.loc['MMseqs2 with 2 Threads, 7.5 Sensitivity',
             'Elapsed (wall clock) time (Seconds)']

str(datetime.timedelta(seconds=seconds32))
str(datetime.timedelta(seconds=seconds2))

