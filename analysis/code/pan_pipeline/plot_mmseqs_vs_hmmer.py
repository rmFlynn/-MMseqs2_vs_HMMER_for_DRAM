"""Download format data and make a plot"""
import sys
import os
import pandas as pd
import numpy as np
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
from multiprocessing import Pool
import seaborn as sns
from shared_tools import OUT_PLOT_PATH

sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from output_tools.plot_tools import LIGHT_COLOR
from vogdb_by_vogdb.shared_tools import \
    OUT_DATA_PATH as VOGDB_BY_VOGDB_OUT_DATA_PATH
from shale_by_vogdb.shared_tools import \
    OUT_DATA_PATH as SHALE_BY_VOGDB_OUT_DATA_PATH
from genome15_by_vogdb.shared_tools import \
    OUT_DATA_PATH as GENOME15_BY_VOGDB_OUT_DATA_PATH
from genome15_by_dbcan.shared_tools import \
    OUT_DATA_PATH as GENOME15_BY_DBCAN_OUT_DATA_PATH
from cazydb_by_dbcan.shared_tools import \
    OUT_DATA_PATH as CAZYDB_BY_DBCAN_OUT_DATA_PATH

data_paths = {
    'vogdb_by_vogdb': VOGDB_BY_VOGDB_OUT_DATA_PATH,
    'shale_by_vogdb': SHALE_BY_VOGDB_OUT_DATA_PATH,
    'genome15_by_vogdb': GENOME15_BY_VOGDB_OUT_DATA_PATH,
    'genome15_by_dbcan': GENOME15_BY_DBCAN_OUT_DATA_PATH,
    'cazydb_by_dbcan': CAZYDB_BY_DBCAN_OUT_DATA_PATH
}

def read_stats(key, path):
    data = pd.read_pickle(os.path.join(
        path, 'mmseqs_vs_hmmer_comp_stats.pkl'))
    data['Pipeline'] = key
    return data

data = pd.concat([read_stats(k, data_paths[k]) for k in data_paths])
# limite to 7.5 for mmseqs_sensitivity
data = data[data['sensitivity'] == 7.5]
# data = data[(data['e-value'] >= 1e-16) & (data['e-value'] <= 1e-14)]
count_cols = ['In MMseqs, Not In HMMER', 'in HMMER, Not In MMseqs', 'In Both']
data['sum'] = data[count_cols].sum(1)
for i in count_cols:
    data[i] = (data[i] / data['sum']) * 100

data = pd.melt(
    data,
    value_vars=count_cols,
    id_vars=[i for i in data.columns if i not in
             count_cols])
data.rename(columns={'variable': "Annotation",
                      'value': "Percent of Total Observations"},
             inplace=True)

data['e-value_exp'] = np.log10(data['e-value'])

pipes = data['Pipeline'].unique()
numplots = len(pipes)
pal = sns.cubehelix_palette(numplots,
                              light=0.7, start=.1,
                              rot=-1.75)
sns.set_style("whitegrid")
plt.clf()
fig, axs = plt.subplots(figsize=(18, 6), ncols=numplots)
for i in range(numplots):
    sns.lineplot(style='Annotation',
                 y="Percent of Total Observations",
                 x='e-value_exp',
                 color=pal[i],
                 ax=axs[i],
                 data=data[data['Pipeline'] == pipes[i]],
                 markers=True
                 )
    axs[i].legend(title='MMseqs Sensitivity Setting')
    axs[i].set_ylim([data['e-value_exp'].min(), 100])
    axs[i].set_title(pipes[i])
    if i == 0:
        axs[0].legend(loc='lower center', ncol=3,
                  bbox_to_anchor=(3, -0.25))
    else:
        axs[i].set_ylabel('')
        axs[i].yaxis.set_ticklabels([])
        axs[i].legend([],[], frameon=False)
    if i != 3:
        axs[i].set_xlabel("MMseqs E-Value Exponent")

fig.subplots_adjust(wspace=0.01, bottom=0.2)
fig.patch.set_facecolor(LIGHT_COLOR)
fig.savefig(os.path.join(
    OUT_PLOT_PATH, 'mmseqs_vs_hmmer_comparison.png'))

def main():
    pass


if __name__ == '__main__':
    main()





