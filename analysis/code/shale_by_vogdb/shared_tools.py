"""All functions shared across this analysis"""

import sys
from io import StringIO
from subprocess import Popen, PIPE
import pandas as pd

HMMER_FILE = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/shale_by_vogdb/"\
    "hmmer/shale_by_vogdb_hmmer_results_2021_06_18_09/hmmer_results.b6"
MMSEQS_SWEEP_OUTPUT = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/shale_by_vogdb/"\
    "mmseqs/shale_by_vogdb_mmseqs_results_2021_06_17_18/sweep_output"
RAW_PROTEINS = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/shale/genes.faa"
OUT_DATA_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/data/"\
    "shale_by_vogdb/"
OUT_LATEX_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/"\
    "tables/shale_by_vogdb"
OUT_PLOT_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/"\
    "plots/shale_by_vogdb"

MMSEQS_ARGS = {'strip_from_annotations':['.msa'],
               'take_only_one':True}
HMMER_ARGS = {'strip_from_annotations':[],
              'file_name': HMMER_FILE,
              'take_only_one':True}


