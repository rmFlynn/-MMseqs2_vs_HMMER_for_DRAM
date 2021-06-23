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

def parse_hmmsearch_domtblout(file_name:str, strip_from_annotations,
                              take_only_one:bool, filter_dbcan_under=False,
                              filter_regex=None, limit_e=1e-15,
                              pct_cover=0.35):
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
    results = results[results['pct_cover'] > pct_cover]
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
    if filter_dbcan_under:
        results['annotation'] = \
            results['annotation'].str.split('_', expand=True)[0]
    if filter_regex is not None:
        results = results[results['annotation'].apply(
            lambda x: re.match(filter_regex, x))]
    return results

def read_mmseqs_results(file_name:str, strip_from_annotations,
                        take_only_one:bool, filter_dbcan_under=False,
                        filter_regex=None, pct_cover=0.35):
    results = pd.read_csv(file_name, sep='\t',
                          header=None, names=BOUTFMT6_COLUMNS)
    # filter percent coverage
    results['pct_cover'] = results.apply(
        lambda x: (x['tEnd'] - x['tStart']) / x['alnLen'],
        axis=1)
    results = results[results['pct_cover'] > pct_cover]
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
    if take_only_one:
        results = results.sort_values('e-value').groupby('ProteinID').\
            first().reset_index()
    # remove everything after under score probably only for dbcan results
    if filter_dbcan_under:
        results['annotation'] = \
            results['annotation'].str.split('_', expand=True)[0]
    # filter based on regex probly only for dbcan results
    if filter_regex is not None:
        results = results[results['annotation'].apply(
            lambda x: re.match(filter_regex, x))]
    return results

def list_proteins(raw_proteins):
    """
    :returns: a count of all proteins
    """
    part1 = Popen(["cat", raw_proteins], stdout=PIPE)
    part2 = Popen(['grep',  '^>'], stdin=part1.stdout, stdout=PIPE)
    all_proteins = pd.read_csv(
        StringIO(
            part2.communicate()[0].decode('utf-8')
        ),
        header=None,
        delimiter=r"\n"
    )
    all_proteins.columns = ['proteins']
    all_proteins['proteins'] = all_proteins['proteins'].str.replace(">", "")
    all_proteins['proteins'] = all_proteins['proteins'].str.\
        split(' ', expand=True)[0]
    return all_proteins['proteins'].unique()

def make_hmmer_vs_mmseqs_df(hmmer, mmseqs):
    hmmer["hmmer_pred"] = 1
    hmmer.rename(columns={"e-value":"hmmer_e-value"}, inplace=True)
    mmseqs["mmseqs_pred"] = 1
    mmseqs.rename(columns={"e-value":"mmseqs_e-value"}, inplace=True)
    data = pd.merge(hmmer, mmseqs, how='outer', on=['ProteinID', 'annotation'])
    data['mmseqs_pred'].fillna(0, inplace=True)
    data['hmmer_pred'].fillna(0, inplace=True)
    return data


def list_annotations(raw_alignments, strip_from_annotations:list):
    """
    :param file_name: The raw aliments file
    :returns: the count of all vogs
    """
    all_anos = pd.read_csv(
        StringIO(
            Popen(
                ["tar", "-tvf", raw_alignments],
                stdout=PIPE).communicate()[0].decode('utf-8')
        ),
        header=None,
    delimiter=r"\s+"
    )
    all_anos =  all_anos[[5]]
    all_anos.columns = ['annotation']

    for ano in strip_from_annotations:
        all_anos['annotation'] = \
            all_anos['annotation'].str.replace(ano, "", regex=False)
    return all_anos['annotation'].unique()


def set_up_comparison(data, comp, on):
    """
    Used in the test_id_validity process

    :param data:
    :param comp:
    :param on:
    :returns:
    """
    data = data.copy()
    comp['is in in'] = 1
    data['is in out'] = 1
    comp_df = pd.merge(data, comp, on=on, how='outer')
    comp_df['is in out'].fillna(0, inplace=True)
    comp_df['is in in'].fillna(0, inplace=True)
    return comp_df


def set_up_comparison(data, comp, on):
    data = data.copy()
    comp['is in in'] = 1
    data['is in out'] = 1
    comp_df = pd.merge(data, comp, on=on, how='outer')
    comp_df['is in out'].fillna(0, inplace=True)
    comp_df['is in in'].fillna(0, inplace=True)
    return comp_df

def set_up_comparison(data, comp, on):
    data = data.copy()
    comp['is in in'] = 1
    data['is in out'] = 1
    comp_df = pd.merge(data, comp, on=on, how='outer')
    comp_df['is in out'].fillna(0, inplace=True)
    comp_df['is in in'].fillna(0, inplace=True)
    return comp_df
