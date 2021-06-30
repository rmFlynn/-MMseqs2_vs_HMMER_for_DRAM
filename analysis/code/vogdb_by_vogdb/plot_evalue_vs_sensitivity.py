"""Download format data and make a plot"""
import sys
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from shared_tools import OUT_PLOT_PATH, OUT_DATA_PATH
sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from output_tools.plot_tools import LIGHT_COLOR

def make_evalues_vs_sensitivty_by_mmseqs_sensitivity_plot(all_stats):
    all_mmseqs = all_stats[all_stats['name'] != 'hmmer_default'].copy()
    all_mmseqs['evalue_exp'] = np.log10(all_mmseqs['evalue'])
    plt.clf()
    fig, ax = plt.subplots(figsize=(7, 7))
    sns.lineplot(hue='mmseqs_sensitivity', y='sensitivity', x='evalue_exp',
                 palette=\
                 sns.cubehelix_palette(light=0.7, start=.1, rot=-.75,
                                       as_cmap=True),
                 ax=ax, data=all_mmseqs)
    ax.set_ylabel("Calculated Sensitivity")
    ax.set_xlabel("MMseqs e-value Setting")
    ax.legend(title='MMseqs Sensitivity Setting')
    ax.grid(False)
    fig.patch.set_facecolor(LIGHT_COLOR)
    fig.savefig(os.path.join(
        OUT_PLOT_PATH, 'evalues_vs_sensitivty_by_mmseqs_sensitivity.png'))

def make_evalues_vs_specificity_by_mmseqs_sensitivity_plot(all_stats):
    all_mmseqs = all_stats[all_stats['name'] != 'hmmer_default'].copy()
    all_mmseqs['evalue_exp'] = np.log10(all_mmseqs['evalue'])
    plt.clf()
    fig, ax = plt.subplots(figsize=(7, 7))
    sns.lineplot(hue='mmseqs_sensitivity', y='specificity', x='evalue_exp',
                 palette=\
                 sns.cubehelix_palette(light=0.7, start=.1, rot=-.75,
                                       as_cmap=True),
                 ax=ax, data=all_mmseqs)
    ax.set_xlabel("MMseqs e-value Setings")
    ax.set_ylabel("Calculated Specificity")
    ax.set_xlabel("MMseqs e-value Setting")
    ax.legend(title='MMseqs Sensitivity Setting')
    ax.grid(False)
    fig.patch.set_facecolor(LIGHT_COLOR)
    fig.savefig(os.path.join(
        OUT_PLOT_PATH, 'evalues_vs_specificity_by_mmseqs_sensitivity.png'))

def make_roc_by_mmseqs_sensitivity_plot(all_stats):
    all_stats['tp_rate'] = all_stats['sensitivity']
    all_stats['fp_rate'] = 100 - all_stats['specificity']
    all_mmseqs = all_stats[all_stats['name'] != 'hmmer_default'].copy()
    all_hmmer = all_stats[all_stats['name'] == 'hmmer_default'].copy()
    plt.clf()
    fig, ax = plt.subplots(figsize=(7, 7))
    sns.lineplot(y='tp_rate', x='fp_rate', hue='mmseqs_sensitivity',
                 palette=\
                 sns.cubehelix_palette(light=0.7, start=.1, rot=-.75,
                                       as_cmap=True),
                 ax=ax, data=all_mmseqs)
    ax.scatter(x=all_hmmer['fp_rate'], y=all_hmmer['tp_rate'], c="orange")
    ax.annotate("HMMER3", xy=(all_hmmer['fp_rate'], all_hmmer['tp_rate']),
                horizontalalignment='right')
    for _, row in all_mmseqs[all_mmseqs['evalue'] == 1e-15].iterrows():
        ax.scatter(x=row['fp_rate'], y=row['tp_rate'], c="blue")
    ax.set_xlabel("FP Rate")
    ax.set_ylabel("TP Rate")
    ax.plot([-200, -200], [0, 0], linewidth=2)
    ax.legend(title='MMseqs Sensitivity Setting')
    ax.grid(False)
    fig.patch.set_facecolor(LIGHT_COLOR)
    fig.savefig(os.path.join(OUT_PLOT_PATH, 'roc_by_mmseqs_sensitivity.png'))

def make_roc_by_mmseqs_evalue_plot(all_stats):
    all_stats['tp_rate'] = all_stats['sensitivity']
    all_stats['fp_rate'] = 100 - all_stats['specificity']
    all_stats['evalue_exp'] = np.log10(all_stats['evalue'])
    plt.clf()
    fig, ax = plt.subplots(figsize=(7, 7))
    all_mmseqs = all_stats[all_stats['name'] != 'hmmer_default'].copy()
    all_hmmer = all_stats[all_stats['name'] == 'hmmer_default'].copy()
    sns.lineplot(y='tp_rate', x='fp_rate', hue='evalue_exp',
                 palette=\
                 sns.cubehelix_palette(light=0.7, start=.1, rot=-.75,
                                       as_cmap=True),
                 ax=ax, data=all_mmseqs)
    ax.scatter(x=all_hmmer['fp_rate'], y=all_hmmer['tp_rate'], c="orange")
    ax.annotate("HMMER3", xy=(all_hmmer['fp_rate'], all_hmmer['tp_rate']),
                horizontalalignment='right')
    # NOTE: 1e-15 is marked in blue
    for _, row in all_mmseqs[all_mmseqs['evalue'] == 1e-15].iterrows():
        ax.scatter(x=row['fp_rate'], y=row['tp_rate'], c="blue")

    ax.set_xlabel("FP Rate")
    ax.set_ylabel("TP Rate")
    ax.legend(title='MMseqs e-value Setting')
    ax.grid(False)
    fig.patch.set_facecolor(LIGHT_COLOR)
    fig.savefig(os.path.join(OUT_PLOT_PATH, 'roc_by_mmseqs_evalue.png'))

def main():
    all_stats = pd.read_pickle(os.path.join(OUT_DATA_PATH,
                                            'gold_standard_comp_stats.pkl'))
    all_stats['evalue'] = all_stats['name'].apply(
        lambda x: float(x.split('_')[1]) if x != 'hmmer_default' else 1e-15)
    all_stats['mmseqs_sensitivity'] = all_stats['name'].apply(
        lambda x: float(x.split('_')[3]) if x != 'hmmer_default' else np.nan)

    sns.set_style("whitegrid")
    make_evalues_vs_sensitivty_by_mmseqs_sensitivity_plot(all_stats)
    make_evalues_vs_specificity_by_mmseqs_sensitivity_plot(all_stats)
    make_roc_by_mmseqs_sensitivity_plot(all_stats)
    make_roc_by_mmseqs_evalue_plot(all_stats)

if __name__ == '__main__':
    main()





