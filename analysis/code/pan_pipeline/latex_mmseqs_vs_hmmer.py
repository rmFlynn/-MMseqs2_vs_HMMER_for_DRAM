"""Download format data and make a plot"""
import sys
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from shared_tools import OUT_LATEX_PATH

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
count_pct_cols = [i + " (Count, Percent)" for i in count_cols]

data.rename(columns={i:j for i, j in zip(count_cols, count_pct_cols)},
            inplace=True)
data.columns
data['sum'] = data[count_pct_cols].sum(1)

for i in count_pct_cols:
    data[i] = data.apply(lambda x: "{}, {:.2f}%".format(x[i], (x[i] / x['sum']) * 100), 1)

data.sort_values('e-value', inplace=True, ascending=False)
data_sets = {i: data[data['e-value'] == i].copy()
             for i in data['e-value'].unique()}

for i in data_sets:
    data_sets[i] = data_sets[i][count_pct_cols + ['Pipeline']]
    data_sets[i].set_index('Pipeline', drop=True, inplace=True)
    data_sets[i] = data_sets[i].T

for i in data_sets:
    with open(os.path.join(
        OUT_LATEX_PATH,
        "evalue_{:.1}_mmseqs_vs_hmmer_s7.5.tex".format(i)), "w") as f:
        f.write(data_sets[i].to_latex(column_format='r|ccccc'))



def main():
    pass


if __name__ == '__main__':
    main()





