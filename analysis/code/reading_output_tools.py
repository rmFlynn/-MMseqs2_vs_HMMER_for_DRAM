"""shared tools for the reading of vog files"""
import sys
from io import StringIO
from multiprocessing import Pool
from subprocess import Popen, PIPE
import pandas as pd
sys.path.insert(0, '/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/src/')
from shared import BOUTFMT6_COLUMNS

HMMSCAN_ALL_COLUMNS = [
    'query_id', 'query_ascession', 'query_length', 'target_id',
    'target_ascession', 'target_length', 'full_evalue', 'full_score',
    'full_bias', 'domain_number', 'domain_count', 'domain_cevalue',
    'domain_ievalue', 'domain_score', 'domain_bias', 'target_start',
    'target_end', 'alignment_start', 'alignment_end', 'query_start',
    'query_end', 'accuracy', 'description'
]

HMMSCAN_COLUMN_TYPES = [str, str, int, str, str, int, float, float, float,
                        int, int, float, float, float, float, int, int, int,
                        int, int, int, float, str]
PCT_COVER = 0.35

def parse_hmmsearch_domtblout_parent(file_name:str, strip_from_annotations,
                                     take_only_one:bool, limit_e=1e-15):
    df_lines = list()
    for line in open(file_name):
        if not line.startswith('#'):
            line = line.split()
            line = line[:22] + [' '.join(line[22:])]
            df_lines.append(line)
    # Name columns
    results = pd.DataFrame(df_lines, columns=HMMSCAN_ALL_COLUMNS)
    # filter percent coverage
    results[["target_start", "target_end", "target_length"]] = \
        results[["target_start", "target_end", "target_length"]].astype(float)
    results['pct_cover'] = results.apply(
        lambda x: (x['target_end'] - x['target_start']) / x['target_length'],
        axis=1)
    results = results[results['pct_cover'] > PCT_COVER]
    # Slim the df
    results = results[['query_id', 'target_id', 'full_evalue']]
    # rename columns
    results.columns = ['ProteinID', 'annotation', 'e-value']
    # Retype columns
    results = results.astype({'ProteinID': str, 'annotation': str, 'e-value':float})

    for stript in strip_from_annotations:
        results['annotation'] = results['annotation'].\
            str.replace(stript, '', regex=False)
    # Take only the best evalues for each
    if take_only_one :
        results = results.sort_values('e-value').groupby('ProteinID').first().\
            reset_index()
    # filter evalue to match defalt
    results = results[results['e-value'] < limit_e]
    return results

def read_mmseqs_results_parent(file_name:str, strip_from_annotations,
                        take_only_one:bool, pct_cover=0.35):
    results = pd.read_csv(file_name, sep='\t',
                          header=None, names=BOUTFMT6_COLUMNS)
    # filter percent coverage
    results['pct_cover'] = results.apply(
        lambda x: (x['tEnd'] - x['tStart']) / x['alnLen'],
        axis=1)
    results = results[results['pct_cover'] > PCT_COVER]
    # Remove file type extension
    for stript in strip_from_annotations:
        results['tId'] = results['tId'].\
            str.replace(stript, '', regex=False)
    # Slim the df
    results = results[['qId', 'tId', 'eVal']]
    results.rename(columns={'qId':'ProteinID',
                             'tId':'annotation',
                             'eVal':'e-value'},
                   inplace=True)
    # Take only the best evalues for each
    if take_only_one :
        results = results.sort_values('e-value').groupby('ProteinID').\
            first().reset_index()
    return results

