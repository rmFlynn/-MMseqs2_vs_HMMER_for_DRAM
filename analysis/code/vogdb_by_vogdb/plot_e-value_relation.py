"""Download format data and make a plot"""
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
from shared_tools import parse_hmmsearch_domtblout, read_mmseqs_results, \
    make_hmmer_vs_mmseqs_df, OUT_PLOT_PATH, LIGHT_COLOR, MMSEQS_SWEEP_OUTPUT

def make_plot_data():
    hmmer = parse_hmmsearch_domtblout(limit_e=1e-2)
    mmseqs = read_mmseqs_results(os.path.join(
        MMSEQS_SWEEP_OUTPUT, 'evalue_0.01_sens_7.500000'))

    mmseqs['e-value'].min()
    mmseqs['e-value'].max()
    hmmer['e-value'].min()
    hmmer['e-value'].max()

    data = make_hmmer_vs_mmseqs_df(hmmer, mmseqs)
    data.dropna(inplace=True)

    data['mmseqs_e-value_exp'] = np.log10(data['mmseqs_e-value'])
    data['hmmer_e-value_exp'] = np.log10(data['hmmer_e-value'])
    return data

def make_scater(data, zoom=50, evalue=1e-15):
    plt.clf()
    fig, ax = plt.subplots(figsize=(7, 7))
    # ax.scatter( data['mmseqs_e-value'], data['hmmer_e-value'])
    from sklearn.linear_model import LinearRegression
    reg = LinearRegression().fit(data[['hmmer_e-value']].values,
                                 data['mmseqs_e-value'].values)
    reg_x = [hmmer['hmmer_e-value'].min(), hmmer['hmmer_e-value'].max()]
    reg_y = reg.predict(np.array([reg_x]).T)
    data['mmseqs_e-value_exp'] = np.log10(data['mmseqs_e-value'])
    data['hmmer_e-value_exp'] = np.log10(data['hmmer_e-value'])
    data['Significant 1e-15'] = data.apply(
        lambda x: \
        "HMMER" \
        if x['hmmer_e-value'] < evalue and x['mmseqs_e-value'] >= evalue else \
        "MMseqs" \
        if x['hmmer_e-value'] >= evalue and x['mmseqs_e-value'] < evalue else \
        'Both'
        if x['hmmer_e-value'] < evalue and x['mmseqs_e-value'] < evalue else \
        'Neither'
        , axis=1)
    # ax.set_xscale('log', base= 10)
    # ax.set_yscale('log', base= 10)
    ax.set_xlim(-zoom, 0)
    ax.set_ylim(-zoom, 0)
    sns.scatterplot(y='mmseqs_e-value_exp', x='hmmer_e-value_exp',
                    alpha=0.05, palette=sns.cubehelix_palette(4, rot=1,
                                                              light=0.65),
                    hue='Significant 1e-15', ax=ax, data=data,
                    legend=True)
    # ax.scatter( data['mmseqs_e-value'], data['hmmer_e-value'])
    ax.set_xlabel("HMMER e-value (log base 10)")
    ax.set_ylabel("MMseqs2 e-value (log base 10)")
    fig.patch.set_facecolor(LIGHT_COLOR)
    fig.savefig(os.path.join(
        OUT_PLOT_PATH, 'evalue_relation_log_%i.png' %zoom))

if __name__ == "__main__":
    data = make_plot_data()
    make_scater(data, zoom=50)
    make_scater(data, zoom=100)

# plt.clf()
# fig, ax = plt.subplots(figsize=(7, 7))
# # ax.scatter( data['mmseqs_e-value'], data['hmmer_e-value'])
# sns.set_palette("Set2_r")
# data['mmseqs_e-value_exp'] = np.log10(data['mmseqs_e-value'])
# data['hmmer_e-value_exp'] = np.log10(data['hmmer_e-value'])
# data = data[data['mmseqs_e-value_exp'] > -100]
# data = data[data['hmmer_e-value_exp'] > -100]
# sns.histplot(y='mmseqs_e-value_exp', x='hmmer_e-value_exp',
#                  ax=ax, data=data,kde=True
#                 )
# # ax.scatter( data['mmseqs_e-value'], data['hmmer_e-value'])
# ax.set_xlabel("HMMER e-value (log base 10)")
# ax.set_ylabel("MMseqs2 e-value (log base 10)")
# fig.patch.set_facecolor(LIGHT_COLOR)
# fig.savefig('plots/evalue_relation_log_kde.png')

