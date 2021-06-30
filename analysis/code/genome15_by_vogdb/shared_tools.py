"""All functions shared across this analysis"""
import sys
from io import StringIO
from subprocess import Popen, PIPE
import pandas as pd

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
DARKER_COLOR = (220/255, 220/255, 220/255)
LIGHT_COLOR = (240/255, 240/255, 240/255)

MMSEQS_ARGS = {'strip_from_annotations':['.msa'],
               'take_only_one':True}
HMMER_ARGS = {'strip_from_annotations':[],
              'file_name': HMMER_FILE,
              'take_only_one':True}


