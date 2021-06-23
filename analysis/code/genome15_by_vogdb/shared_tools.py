"""All functions shared across this analysis"""
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
    "genome15_by_vogdb/hmmer/genome15_by_vogdb_hmmer_results_2021_06_18_08/hmmer_results.b6"
MMSEQS_SWEEP_OUTPUT = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/"\
    "genome15_by_vogdb/mmseqs/genome15_by_vogdb_mmseqs_results_2021_06_17_18/sweep_output"
RAW_PROTEINS = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/genome15/"\
    "genes.faa"
OUT_DATA_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/data/"\
    "genome15_by_vogdb/"
OUT_LATEX_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/"\
    "tables/genome15_by_vogdb"
OUT_PLOT_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/"\
    "plots/genome15_by_vogdb"

MMSEQS_ARGS = {'strip_from_annotations':['.msa'],
               'take_only_one':True}
HMMER_ARGS = {'strip_from_annotations':[],
              'file_name': HMMER_FILE,
              'take_only_one':True}

def parse_hmmsearch_domtblout():
    return parse_hmmsearch_domtblout_parent(
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
        delimiter=r"\n",
        engine='python'
    )
    all_proteins.columns = ['proteins']
    all_proteins['proteins'] = all_proteins['proteins'].str.replace(">", "")
    all_proteins['proteins'] = all_proteins['proteins'].str.\
        split(' ', expand=True)[0]
    return all_proteins['proteins'].unique()
