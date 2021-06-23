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

BOUTFMT6_COLUMNS_MMSEQS = ['ProteinID', 'dbcan_id', 'seqIdentity', 'alnLen',
                           'mismatchCnt', 'gapOpenCnt', 'qStart', 'qEnd',
                           'tStart', 'tEnd', 'e-value', 'bitScore']
HMMSCAN_ALL_COLUMNS = [
    'query_id', 'query_ascession', 'query_length', 'target_id',
    'target_ascession', 'target_length', 'full_evalue', 'full_score',
    'full_bias', 'domain_number', 'domain_count', 'domain_cevalue',
    'domain_ievalue', 'domain_score', 'domain_bias', 'target_start',
    'target_end', 'alignment_start', 'alignment_end', 'query_start',
    'query_end', 'accuracy', 'description']
HMMSCAN_COLUMN_TYPES = [str, str, int, str, str, int, float, float, float,
                        int, int, float, float, float, float, int, int, int,
                        int, int, int, float, str]

def parse_hmmsearch_domtblout(file_name:str, pct_cover=0.35):
    df_lines = list()
    for line in open(file_name):
        if not line.startswith('#'):
            line = line.split()
            line = line[:22] + [' '.join(line[22:])]
            df_lines.append(line)
    results = pd.DataFrame(df_lines, columns=HMMSCAN_ALL_COLUMNS)
    results[["target_start", "target_end", "target_length"]] = \
        results[["target_start", "target_end", "target_length"]].astype(float)
    results['pct_cover'] = results.apply(
        lambda x: (x['target_end'] - x['target_start']) / x['target_length'],
        axis=1)
    # filter percent coverage
    results = results[results['pct_cover'] > pct_cover]
    results = results[['query_id', 'target_id', 'full_evalue']]
    results.columns = ['ProteinID', 'dbcan_id', 'e-value']
    results = results.astype({'ProteinID': str, 'dbcan_id': str, 'e-value':float})
    results['dbcan_id'] = results['dbcan_id'].\
        str.replace('.hmm', '', regex=False)
    # filter evalue to match curent
    results = results[results['e-value'] < 1e-15]
    # # Take only the best
    # results = results.sort_values('e-value').groupby('ProteinID').first().\
    #     reset_index()
    return results

def read_mmseqs_results(path, pct_cover=0.35):
    results = pd.read_csv(path, sep='\t',
                          header=None, names=BOUTFMT6_COLUMNS_MMSEQS)
    results['pct_cover'] = results.apply(
        lambda x: (x['tEnd'] - x['tStart']) / x['alnLen'],
        axis=1)
    results = results[results['pct_cover'] > pct_cover]
    results['dbcan_id'] = results['dbcan_id'].\
        str.replace('.aln', '', regex=False)
    results = results[['ProteinID', 'dbcan_id', 'e-value']]
    # # Take only the best
    # results = results.sort_values('e-value').groupby('ProteinID').first().\
    #     reset_index()
    return results

def make_combo_df(hmmer, mmseqs):
    hmmer["hmmer_pred"] = 1
    hmmer.rename(columns={"e-value":"hmmer_e-value"}, inplace=True)
    mmseqs["mmseqs_pred"] = 1
    mmseqs.rename(columns={"e-value":"mmseqs_e-value"}, inplace=True)
    data = pd.merge(hmmer, mmseqs, how='outer', on=['ProteinID', 'dbcan_id'])
    data['mmseqs_pred'].fillna(0, inplace=True)
    data['hmmer_pred'].fillna(0, inplace=True)
    return data

def make_dif_df(data, name):
    # this is a confusion_matrix but hmmer is true
    differ = confusion_matrix(data['hmmer_pred'], data['mmseqs_pred'])
    differ_df = {}
    _, differ_df['In MMseqs, Not In HMMER'], \
        differ_df['in HMMER, Not In MMseqs'], differ_df['both'] = \
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


# hmmer_raw = \
#     parse_hmmsearch_domtblout("../hmmer/hmmer_cazydb_results.b6")
# hmmer_raw.to_pickle("hmmer_raw.pkl")
# hmmer == hmmer_raw.copy()
hmmer = pd.read_pickle("hmmer_raw.pkl")
hmmer['True ID'] = hmmer['ProteinID'].str.split('|')
hmmer.rename(columns={'dbcan_id': 'Predicted ID'}, inplace=True)
hmmer['Predicted ID'] = hmmer['Predicted ID'].str.split('_', expand=True)[0]
hmmer['True ID'] = hmmer['True ID'].\
    apply(lambda x: [i.split('_')[0] for i in x])
hmmer['True Pred'] = hmmer.apply(
    lambda x: 1 if x['Predicted ID'] in x['True ID'] else 0, axis=1)
hmmer = hmmer.sort_values('e-value').\
    groupby(['ProteinID', 'Predicted ID']).\
    first().\
    reset_index()
hmmer.rename(columns={'e-value': 'HMMER E-Value'}, inplace=True)
hmmer['HMMER Pred'] = 1
hmmer[hmmer['True Pred'] == 0][['Predicted ID', 'True ID']]


# hmmer['True ID'] = hmmer['ProteinID'].str.split('|', expand=True)[1]
# hmmer.rename(columns={'dbcan_id': 'Predicted ID'}, inplace=True)
# # NOTE that the ids were simplified
# hmmer['True ID'] = hmmer['True ID'].apply(lambda x: str(x).split('_')[0])
# # NOTE dropping a small number of obs that don't conform to ProteinID format
# #      used to extract ids
# hmmer = hmmer[hmmer['True ID'] != 'None']
# hmmer['Predicted ID'] = hmmer['Predicted ID'].apply(lambda x: x.split('_')[0])
# hmmer = hmmer.sort_values('e-value').\
#     groupby(['Predicted ID', 'True ID']).\
#     first().\
#     reset_index()
#
# hmmer[(hmmer['True ID'].str.match(r'[A-z]+\d+')) == False]
# # NOTE many of these are from other data sets they will be treated as wrong
# cant_match = hmmer[(hmmer['Predicted ID'].\
#     str.\
#     match(r'[A-z]+\d+')) == False][['Predicted ID', 'True ID']].\
#     drop_duplicates()
# cant_match['Predicted ID'].unique()
# hmmer.rename(columns={'e-value': 'HMMER E-Value'}, inplace=True)
# hmmer['HMMER Pred'] = 1
#
# hmmer[hmmer['Predicted ID'] != \
#       hmmer['True ID']][['Predicted ID', 'True ID']]

# mmseqs_dir = '../mmseqs/cazydb_results_2021_06_15_14/sweep_output/'
# mmseqs_raw = read_mmseqs_results(
#     os.path.join(mmseqs_dir, 'evalue_1e-15_sens_7.500000'))
# mmseqs_raw.to_pickle("mmseqs_raw.pkl")
# mmseqs == mmseqs_raw.copy()
mmseqs = pd.read_pickle("mmseqs_raw.pkl")
mmseqs['True ID'] = mmseqs['ProteinID'].str.split('|')
mmseqs.rename(columns={'dbcan_id': 'Predicted ID'}, inplace=True)
mmseqs['Predicted ID'] = mmseqs['Predicted ID'].str.split('_', expand=True)[0]
mmseqs['True ID'] = mmseqs['True ID'].\
    apply(lambda x: [i.split('_')[0] for i in x])
mmseqs['True Pred'] = mmseqs.apply(
    lambda x: 1 if x['Predicted ID'] in x['True ID'] else 0, axis=1)
mmseqs = mmseqs.sort_values('e-value').\
    groupby(['ProteinID', 'Predicted ID']).\
    first().\
    reset_index()
mmseqs.rename(columns={'e-value': 'MMseqs E-Value'}, inplace=True)
mmseqs['MMseqs Pred'] = 1
mmseqs[mmseqs['True Pred'] == 0][['Predicted ID', 'True ID']]
# true_ids = mmseqs['ProteinID'].\
#     apply(lambda x: '|'.join(x.split('|')[1:-1])).\
#     str.split('|', expand = True)
# mmseqs = pd.melt(
#     pd.concat([ mmseqs, true_ids ], axis=1),
#     value_vars=true_ids.columns,
#     id_vars=mmseqs.columns).\
#     drop('variable', axis=1).\
#     rename(columns={'value':'True ID'}).\
#     reset_index(drop=True)
#
# mmseqs = mmseqs[~mmseqs['True ID'].isnull()]
# # NOTE all this filtering
# mmseqs = mmseqs[mmseqs['True ID'] != '']
# mmseqs = mmseqs[~mmseqs['True ID'].str.match(r'\d+\.\d+\.\d+\.\d+')]
# mmseqs = mmseqs[~mmseqs['True ID'].str.match(r'\d+\.\d+\.\d+\.\-')]
# mmseqs = mmseqs[~mmseqs['True ID'].str.match(r'\d+\.\d+\.\d+\.n1')]
# mmseqs[(mmseqs['True ID'].str.match(r'[A-z]+\d+')) == False]['True ID']
# mmseqs['True ID'] = mmseqs['True ID'].apply(lambda x: x.split('_')[0])
# mmseqs['Predicted ID'] = mmseqs['Predicted ID'].\
#     apply(lambda x: x.split('_')[0])
# for i in mmseqs['ProteinID']:
#     print(i)
#
# mmseqs = mmseqs.sort_values('e-value').\
#     groupby(['Predicted ID', 'True ID']).\
#     first().\
#     reset_index()
# mmseqs[(mmseqs['True ID'].str.match(r'[A-z]+\d+')) == False]['True ID']
# mmseqs[(mmseqs['Predicted ID'].str.\
#         match(r'[A-z]+\d+')) == False]['Predicted ID']

mmseqs.drop('True ID', axis=1, inplace=True)
hmmer.drop('True ID', axis=1, inplace=True)
data = pd.merge(mmseqs, hmmer, on=['ProteinID', 'Predicted ID', 'True Pred'],
                how='outer')
# data['True'] = (data['Predicted ID'] == data['True ID']).astype(int)
data['MMseqs Pred'].fillna(0, inplace=True)
data['HMMER Pred'].fillna(0, inplace=True)

cfmx_mmseqs = confusion_matrix(data['True Pred'], data['MMseqs Pred'])
cfmx_hmmer = confusion_matrix(data['True Pred'], data['HMMER Pred'])
confusion_matrix(data['HMMER Pred'], data['MMseqs Pred'])

cfmx = confusion_matrix(mmseqs['True ID'], mmseqs['Predic ID'])

