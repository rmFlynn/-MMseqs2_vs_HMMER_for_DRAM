"""All functions and objects shared across projects anotated by vogdb"""
import sys
from io import StringIO
from subprocess import Popen, PIPE
import pandas as pd

RAW_ALIGNMENTS = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/"\
    "vog.raw_algs.20210525.tar.gz"

def make_hmmmer_vs_vog_truth_df(hmmer_results, vog_truth):
    hmmer = hmmer_results.copy()
    hmmer["hmmer_pred"] = 1
    vog_truth["true_pred"] = 1
    hmmer.rename(columns={"e-value":"hmmer_e-value"}, inplace=True)
    data = pd.merge(hmmer, vog_truth, how='outer',
                        on=['ProteinID', 'annotation'])
    data['hmmer_e-value'].fillna(10, inplace=True)
    data['hmmer_pred'].fillna(0, inplace=True)
    data['true_pred'].fillna(0, inplace=True)
    return data

def make_mmseqs_vs_vog_truth_df(mmseqs_results, vog_truth):
    mmseqs = mmseqs_results.copy()
    mmseqs["mmseqs_pred"] = 1
    vog_truth["true_pred"] = 1
    mmseqs.rename(columns={"e-value":"mmseqs_e-value"}, inplace=True)
    data = pd.merge(mmseqs, vog_truth, how='outer',
                        on=['ProteinID', 'annotation'])
    data['mmseqs_e-value'].fillna(10, inplace=True)
    data['mmseqs_pred'].fillna(0, inplace=True)
    data['true_pred'].fillna(0, inplace=True)
    return data


def read_vog_truth(file_name):
    data = pd.read_csv(file_name, sep='\t')
    data = pd.concat(
        [data[data.columns[:-1]],
         data['ProteinIDs'].str.split(',', expand=True)],
        axis=1).\
        melt( id_vars=data.columns[:-1],
             var_name="drop",
             value_name='ProteinID').\
        drop(['drop'], axis=1).\
        dropna(subset=['ProteinID'])
    # data.drop_duplicates(subset=['ProteinID'], inplace=True)
    data.rename(columns={'#GroupName':'annotation'}, inplace=True)
    data = data[['ProteinID', 'annotation']]
    return data
