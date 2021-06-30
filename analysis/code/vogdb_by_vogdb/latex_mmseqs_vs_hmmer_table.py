"""Download format data and make a plot"""
import os
import pandas as pd
import numpy as np
from shared_tools import OUT_DATA_PATH, OUT_LATEX_PATH

def make_select_tex(data_path, out_latex_path):
    stats = pd.read_pickle(os.path.join(
        data_path, 'mmseqs_vs_hmmer_comp_stats.pkl'))
    stats = stats[stats['sensitivity'] == 7.5].sort_values(
        'e-value', ascending=False)
    stats = stats[['e-value', 'In MMseqs, Not In HMMER',
           'in HMMER, Not In MMseqs', 'In Both']]
    stats['e-value'] = stats['e-value'].apply(lambda x: "{:.0e}".format(x))

    with open(os.path.join(
        out_latex_path, "mmseqs_vs_hmmer_s7.5.tex"), "w") as f:
        f.write(stats.to_latex(index=False, column_format='r|lll'))

if __name__ == '__main__':
    make_select_tex(OUT_DATA_PATH, OUT_LATEX_PATH)


