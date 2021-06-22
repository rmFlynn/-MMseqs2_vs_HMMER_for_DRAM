"""All functions and objects shared across data bases anotated by vogdb"""
import sys
from io import StringIO
from subprocess import Popen, PIPE
import pandas as pd

RAW_ALIGNMENTS = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/"\
    "vog.raw_algs.20210525.tar.gz"

def list_vogs():
    """
    :param file_name: The raw aliments file
    :returns: the count of all vogs
    """
    all_vogs = pd.read_csv(
        StringIO(
            Popen(
                ["tar", "-tvf", RAW_ALIGNMENTS],
                stdout=PIPE).communicate()[0].decode('utf-8')
        ),
        header=None,
    delimiter=r"\s+"
    )
    all_vogs =  all_vogs[[5]]
    all_vogs.columns = ['annotation']
    all_vogs['annotation'] = all_vogs['annotation'].str.replace(".msa", "", regex=False)
    return all_vogs['annotation'].unique()

def make_hmmer_vs_mmseqs_df(hmmer, mmseqs):
    hmmer["hmmer_pred"] = 1
    hmmer.rename(columns={"e-value":"hmmer_e-value"}, inplace=True)
    mmseqs["mmseqs_pred"] = 1
    mmseqs.rename(columns={"e-value":"mmseqs_e-value"}, inplace=True)
    data = pd.merge(hmmer, mmseqs, how='outer', on=['ProteinID', 'annotation'])
    data['mmseqs_pred'].fillna(0, inplace=True)
    data['hmmer_pred'].fillna(0, inplace=True)
    return data
