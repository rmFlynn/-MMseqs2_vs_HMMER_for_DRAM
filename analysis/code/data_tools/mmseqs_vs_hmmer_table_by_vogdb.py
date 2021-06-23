
import sys
import os
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix
sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from reading_output_tools import parse_hmmsearch_domtblout, \
    read_mmseqs_results
from vog_tools import make_hmmer_vs_mmseqs_df

def make_dif_df(data, name):
    # this is a confusion_matrix but hmmer is true
    differ = confusion_matrix(data['hmmer_pred'], data['mmseqs_pred'])
    differ_df = {}
    _, differ_df['In MMseqs, Not In HMMER'], \
        differ_df['in HMMER, Not In MMseqs'], differ_df['In Both'] = \
        differ.ravel()
    differ_df = pd.DataFrame(
        differ_df, index=[name])
    return differ_df

def split_out_eval_sens(data):
    data['e-value'] = np.vectorize(
        lambda x:float(x.split('_')[1]))(data.index.values)
    data['sensitivity'] = np.vectorize(
        lambda x:float(x.split('_')[3]))(data.index.values)
    return data

def sweep_the_mmseqs(name, path, hmmer, mmseqs_args):
    mmseqs = read_mmseqs_results(os.path.join(path, name), **mmseqs_args)
    return make_dif_df(make_hmmer_vs_mmseqs_df(hmmer, mmseqs), name)

def make_count_table(mmseqs_sweep_output, out_data_path, hmmer_args, mmseqs_args):
    hmmer_results = parse_hmmsearch_domtblout(**hmmer_args)
    count_table = pd.concat([
        sweep_the_mmseqs(name, mmseqs_sweep_output, hmmer_results, mmseqs_args)
        for name in os.listdir(mmseqs_sweep_output)])
    count_table = split_out_eval_sens(count_table)
    count_table.to_csv(os.path.join(out_data_path, 'mmseqs_vs_hmmer_comp_stats.csv'))
    count_table.to_pickle(os.path.join(out_data_path, 'mmseqs_vs_hmmer_comp_stats.pkl'))
