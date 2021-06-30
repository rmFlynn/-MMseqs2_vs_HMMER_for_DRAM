"""Make a plot showing the relationship of the hmmer and mmseqs e-values around the decision point"""

import sys
import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression

sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from reading_output_tools import parse_hmmsearch_domtblout, \
    read_mmseqs_results, make_hmmer_vs_mmseqs_df
from output_tools.plot_tools import LIGHT_COLOR


def make_plot_data(mmseqs_sweep_output, mmseqs_args, hmmer_args):
    hmmer_args['limit_e'] = 1e-2
    hmmer = parse_hmmsearch_domtblout(**hmmer_args)
    mmseqs = read_mmseqs_results(
        os.path.join(mmseqs_sweep_output, 'evalue_0.01_sens_7.500000'),
        **mmseqs_args)

    mmseqs['e-value'].min()
    mmseqs['e-value'].max()
    hmmer['e-value'].min()
    hmmer['e-value'].max()

    data = make_hmmer_vs_mmseqs_df(hmmer, mmseqs)
    data.dropna(inplace=True)

    data['mmseqs_e-value_exp'] = np.log10(data['mmseqs_e-value'])
    data['hmmer_e-value_exp'] = np.log10(data['hmmer_e-value'])
    return data

def make_scater(data, out_plot_path, zoom=50, evalue=1e-15):
    plt.clf()
    fig, ax = plt.subplots(figsize=(7, 7))
    from sklearn.linear_model import LinearRegression
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
    ax.axline([-200, -200], [0, 0], linewidth=2, color='red')
    lin_reg = LinearRegression().fit(
        (lambda d, x, y: d[(-2 >= d[x]) & (d[x] >= -200) & \
                              (-2 >= d[y]) & (d[y] >= -200)][x].\
                             values.reshape(-1, 1))(
                data, "mmseqs_e-value_exp", "hmmer_e-value_exp"),
            (lambda d, x, y: d[(-2 >= d[x]) & (d[x] >= -200) & \
                              (-2 >= d[y]) & (d[y] >= -200)][y].\
                             values.reshape(-1, 1))(
                data, "mmseqs_e-value_exp", "hmmer_e-value_exp")
      )
    lin_reg_cof = lin_reg.coef_[0][0]
    lin_reg_intercept = lin_reg.intercept_[0]
    ax.axline([-200 * lin_reg_cof + lin_reg_intercept, -200],
              [lin_reg_intercept, 0], linewidth=2, color='blue')
    fig.patch.set_facecolor(LIGHT_COLOR)
    fig.savefig(os.path.join(
        out_plot_path, 'evalue_relation_log_%i.png' %zoom))

