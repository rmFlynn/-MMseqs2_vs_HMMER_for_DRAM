"""All functions shared across this analysis"""

import sys
from io import StringIO
from subprocess import Popen, PIPE
import pandas as pd



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
OUT_LATEX_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/"\
    "tables/vogdb_by_vogdb"
OUT_PLOT_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/"\
    "plots/vogdb_by_vogdb"


MMSEQS_ARGS = {'strip_from_annotations':['.msa'],
               'take_only_one':True}
HMMER_ARGS = {'strip_from_annotations':[],
              'file_name': HMMER_FILE,
              'take_only_one':True}


