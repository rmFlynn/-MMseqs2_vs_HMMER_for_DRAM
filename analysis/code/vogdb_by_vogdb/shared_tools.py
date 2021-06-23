"""All functions shared across vog db"""

import sys
from io import StringIO
from subprocess import Popen, PIPE
import pandas as pd

sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from reading_output_tools import parse_hmmsearch_domtblout_parent, \
    read_mmseqs_results_parent
from vog_tools import RAW_ALIGNMENTS, list_vogs, make_hmmer_vs_mmseqs_df


HMMER_FILE = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/"\
    "vogdb_by_vogdb/hmmer/vogdb_by_vogdb_hmmer_results_2021_06_18_09/"\
    "hmmer_results.b6"
MMSEQS_SWEEP_OUTPUT = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/vogdb_by_vogdb/"\
    "mmseqs/vogdb_by_vogdb_mmseqs_t32_results_2021_06_17_18/sweep_output/"
RAW_PROTEINS = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/"\
    "vog.proteins.all.20210525.fa"
VOG_TRUTH_FILE = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/"\
    "vog.members.20210525.tsv.gz"
OUT_DATA_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/data/"\
    "vogdb_by_vogdb/"


MMSEQS_ARGS = {'strip_from_annotations':['.msa'],
               'take_only_one':True}
HMMER_ARGS = {'strip_from_annotations':[],
              'file_name': HMMER_FILE,
              'take_only_one':True}

def parse_hmmsearch_domtblout(limit_e=1e-15):
    return parse_hmmsearch_domtblout_parent(
        limit_e=limit_e,
        **HMMER_ARGS
    )

def read_mmseqs_results(file_name:str):
    return read_mmseqs_results_parent(
        file_name=file_name,
        **MMSEQS_ARGS
    )

def list_vog_protines():
    """
    :returns: a count of all proteins
    """
    part1 = Popen(["cat", RAW_PROTEINS], stdout=PIPE)
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


def read_vog_truth(file_name=VOG_TRUTH_FILE):
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

# Plot and table tools
OUT_LATEX_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/"\
    "tables/vogdb_by_vogdb"
OUT_PLOT_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/"\
    "plots/vogdb_by_vogdb"
DARKER_COLOR = (220/255, 220/255, 220/255)
LIGHT_COLOR = (240/255, 240/255, 240/255)
