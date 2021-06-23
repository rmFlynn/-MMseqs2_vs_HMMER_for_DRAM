"""All functions shared across vog db"""

import sys
import os
import re
from io import StringIO
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np

sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from reading_output_tools import parse_hmmsearch_domtblout_parent, \
    read_mmseqs_results_parent
# from vog_tools import RAW_ALIGNMENTS, list_vogs, make_hmmer_vs_mmseqs_df


HMMER_FILE = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/cazydb_by_dbcan/hmmer/cazydb_by_dbcan_hmmer_results_2021_06_18_10/hmmer_results.b6"
MMSEQS_SWEEP_OUTPUT = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/cazydb_by_dbcan/mmseqs/cazydb_by_dbcan_mmseqs_results_2021_06_18_10/sweep_output/"
RAW_PROTEINS = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/dbcan/CAZyDB.07312020.fa"
RAW_ALIGNMENTS= \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/dbcan/dbCAN-fam-aln-V9.tar.gz"
OUT_DATA_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/data/cazydb_by_dbcan/"


PROTEIN_RE = r"^[A-Z]+\d+$"
MMSEQS_ARGS = {'strip_from_annotations':['.aln'],
               'take_only_one':False}
HMMER_ARGS = {'strip_from_annotations':['.hmm'],
              'file_name': HMMER_FILE,
              'take_only_one':False,
              'filter_dbcan_under':True,# removed in testing
              'filter_regex'=PROTEIN_RE # removed in testing
              }


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


def list_proteins():
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
        delimiter=r"\n",
        engine='python'
    )
    all_proteins.columns = ['proteins']
    all_proteins['proteins'] = all_proteins['proteins'].str.replace(">", "")
    return np.unique(all_proteins['proteins'].values)

def flatten_protein_ids(all_proteins):
    all_proteins = pd.Series(all_proteins).str.\
        split('|', expand = True)
    # NOTE: I checked the first column dose not match the regex
    all_proteins = all_proteins.loc[:, 1:].values.flatten()
    all_proteins = all_proteins[all_proteins != None]
    all_proteins = all_proteins[all_proteins != '']
    return all_proteins

def filter_protein_ids(proteins):
    proteins = np.array([str(x).split('_')[0] for x in proteins])
    keep_proteins = np.array([x for x in  proteins
                              if re.match(PROTEIN_RE, x)])
    drop_proteins = np.array([x for x in  proteins
                              if not re.match(PROTEIN_RE, x)])
    return np.unique(keep_proteins), np.unique(drop_proteins)

def split_out_protein_list(proteins):
    proteins = mmseqs_df['ProteinID']
    protein_df = pd.DataFrame(
        proteins.apply(
            lambda x: filter_protine_ids(x.split('|')[1:])).to_list(),
        index=proteins.index, columns=['keep', 'drop'])
    drop = np.unique(np.concatenate(protein_df['drop'].values))
    return protein_df['keep'], drop

def make_mmseqs_vs_cazydb_truth(df):
    data = df.copy()
    data['annotation'], _ = filter_protine_ids(data['annotation'])
    data['true IDs'], _= split_out_protine_list(data['ProteinID'])
    data['true'] = data.apply(
        lambda x: 1 if x['annotation'] in x['true IDs'] else 0, 1)
    return data

def count_proteins():
    return len(filter_protine_ids(flatten_protein_ids(list_proteins()))[0])

def list_anotations(strip_from_annotations:list):
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

    for ano in strip_from_annotations:
        all_vogs['annotation'] = \
            all_vogs['annotation'].str.replace(ano, "", regex=False)
    return all_vogs['annotation'].unique()


# Plot and table tools
OUT_LATEX_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/tables/cazydb_by_dbcan/"
OUT_PLOT_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/plots/cazydb_by_dbcan"
DARKER_COLOR = (220/255, 220/255, 220/255)
LIGHT_COLOR = (240/255, 240/255, 240/255)
