"""Download format data and make a plot"""
import os
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix
from shared_tools import OUT_DATA_PATH, OUT_LATEX_PATH

def make_latex_stat_rep(df, f):
    f.write(df[['name', 'specificity', 'sensitivity', 'balanced_accuracy']].\
            to_latex(index_names=False, header=False))

def make_latex_cfmx_rep(df, f):
    # f.write('confusion_matrix:\n\n')
    f.write(pd.DataFrame(df['confusion_matrix'],
            columns=['0','1'], index=['0','1']).to_latex())

def make_select_tex():
    all_stats = pd.read_pickle(os.path.join(
        OUT_DATA_PATH, 'gold_standard_comp_stats.pkl'))
    with open(os.path.join(
        OUT_LATEX_PATH, "best_balanced_accuracy_stats.tex"), "w") as f:
        make_latex_stat_rep(
            all_stats.sort_values('balanced_accuracy').iloc[-1], f)

    with open(os.path.join(
        OUT_LATEX_PATH, "best_balanced_accuracy_cfmx.tex"), "w") as f:
        make_latex_cfmx_rep(
            all_stats.sort_values('balanced_accuracy').iloc[-1], f)

    with open(os.path.join(
        OUT_LATEX_PATH, "best_sensitivity_stats.tex"), "w") as f:
        make_latex_stat_rep(
            all_stats.sort_values('sensitivity').iloc[-1], f)

    with open(os.path.join(
        OUT_LATEX_PATH, "best_sensitivity_cfmx.tex"), "w") as f:
        make_latex_cfmx_rep(
            all_stats.sort_values('sensitivity').iloc[-1], f)

    with open(os.path.join(
        OUT_LATEX_PATH, "best_specificity_stats.tex"), "w") as f:
        make_latex_stat_rep(
            all_stats.sort_values('specificity').iloc[-1], f)

    with open(os.path.join(
        OUT_LATEX_PATH, "best_specificity_cfmx.tex"), "w") as f:
        make_latex_cfmx_rep(
            all_stats.sort_values('specificity').iloc[-1], f)

    with open(os.path.join(
        OUT_LATEX_PATH, "hmmer_default_stats.tex"), "w") as f:
        make_latex_stat_rep(
            all_stats[all_stats['name'] == 'hmmer_default'].iloc[0], f)

    with open(os.path.join(
        OUT_LATEX_PATH, "hmmer_default_cfmx.tex"), "w") as f:
        make_latex_cfmx_rep(
            all_stats[all_stats['name'] == 'hmmer_default'].iloc[0], f)

    all_stats['fp'] = all_stats['confusion_matrix'].apply(lambda x: x[0,1])
    all_stats['tp'] = all_stats['confusion_matrix'].apply(lambda x: x[1,1])
    hmmer_cfmx = all_stats[all_stats['name'] == 'hmmer_default']\
        ['confusion_matrix'].values[0]
    all_stats['dif'] = all_stats['confusion_matrix'].\
        apply(lambda x: np.sum(np.abs(x - hmmer_cfmx)))

    all_stats['tp'] = all_stats['confusion_matrix'].apply(lambda x: x[1,1])

    most_similar = all_stats.sort_values('dif').iloc[1]

    best_tp_better_fp = all_stats[all_stats['fp'] < \
              all_stats[all_stats['name'] == 'hmmer_default']['fp'].values[0]]\
        .sort_values('tp').iloc[-1]

    with open(os.path.join(
        OUT_LATEX_PATH, "most_similar_stats.tex"), "w") as f:
        make_latex_stat_rep(most_similar, f)

    with open(os.path.join(
        OUT_LATEX_PATH, "most_similar_cfmx.tex"), "w") as f:
        make_latex_cfmx_rep(most_similar, f)

    with open(os.path.join(
        OUT_LATEX_PATH, "best_tp_less_fp_stats.tex"), "w") as f:
        make_latex_stat_rep(best_tp_better_fp, f)

    with open(os.path.join(
        OUT_LATEX_PATH, "best_tp_less_fp_cfmx.tex"), "w") as f:
        make_latex_cfmx_rep(best_tp_better_fp, f)

def make_all_tex():
    all_stats = pd.read_pickle(os.path.join(
        OUT_DATA_PATH, 'gold_standard_comp_stats.pkl'))
    with open(os.path.join(
        OUT_LATEX_PATH, "All.tex"), "w") as f:
        f.write("\\begin{tabular}{r|l}\n")
        f.write("Confusion Matrix & Description")
        f.write(" \\\\\n")
        for _, j in all_stats.iterrows():
            make_latex_cfmx_rep(j, f)
            f.write("&")
            make_latex_stat_rep(j, f)
            f.write(" \\\\\n")
        f.write("\\end{tabular}")

if __name__ == '__main__':
    make_select_tex()
    make_all_tex()


