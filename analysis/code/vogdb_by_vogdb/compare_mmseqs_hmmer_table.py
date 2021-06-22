"""Download format data and make a plot"""
import sys
import os
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE
from io import StringIO
from matplotlib import pyplot as plt
from sklearn.preprocessing import label_binarize
from sklearn import svm, datasets
from sklearn.metrics import confusion_matrix
from time import sleep
from sklearn.preprocessing import MinMaxScaler
from decimal import Decimal
from multiprocessing import Pool
from shared_tools import parse_hmmsearch_domtblout,\
    read_mmseqs_results, make_hmmer_vs_mmseqs_df, MMSEQS_SWEEP_OUTPUT

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

def sweep_the_mmseqs(name, path, hmmer):
    mmseqs = read_mmseqs_results(os.path.join(path, name))
    return make_dif_df(make_hmmer_vs_mmseqs_df(hmmer, mmseqs), name)

def make_count_table():
    hmmer_results = parse_hmmsearch_domtblout()
    count_table = pd.concat([
        sweep_the_mmseqs(name, MMSEQS_SWEEP_OUTPUT, hmmer_results)
        for name in os.listdir(MMSEQS_SWEEP_OUTPUT)])
    count_table = split_out_eval_sens(genome15_count)
    count_table.to_csv(os.path.join(OUT_DATA_PATH, 'gold_standard_comp_stats.csv'))
    count_table.to_pickle(os.path.join(OUT_DATA_PATH, 'gold_standard_comp_stats.pkl'))

