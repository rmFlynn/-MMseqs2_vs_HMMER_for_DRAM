"""All functions shared across this analysis"""

import sys
import os
import re
from io import StringIO
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np

sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from reading_output_tools import parse_hmmsearch_domtblout, \
    read_mmseqs_results, list_proteins
from dbcan_tools import PROTEIN_RE


HMMER_FILE = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/cazydb_by_dbcan/hmmer/cazydb_by_dbcan_hmmer_results_2021_06_18_10/hmmer_results.b6"
MMSEQS_SWEEP_OUTPUT = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/cazydb_by_dbcan/mmseqs/cazydb_by_dbcan_mmseqs_results_2021_06_18_10/sweep_output/"
RAW_PROTEINS = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/dbcan/CAZyDB.07312020.fa"
OUT_DATA_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/data/cazydb_by_dbcan/"


MMSEQS_ARGS = {'strip_from_annotations':['.aln'],
               'take_only_one':False,
              'filter_dbcan_under':True,# removed in testing
              'filter_regex':PROTEIN_RE # removed in testing
               }
HMMER_ARGS = {'strip_from_annotations':['.hmm'],
              'file_name': HMMER_FILE,
              'take_only_one':False,
              'filter_dbcan_under':True,# removed in testing
              'filter_regex':PROTEIN_RE # removed in testing
              }
MMSEQS_ARGS_TESTING = {'strip_from_annotations':['.aln'],
               'take_only_one':False,
               }
HMMER_ARGS_TESTING = \
    {'strip_from_annotations':['.hmm'],
     'file_name': HMMER_FILE,
     'take_only_one':False,
     }

# Plot and table tools
OUT_LATEX_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/tables/cazydb_by_dbcan/"
OUT_PLOT_PATH = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/analysis/output/plots/cazydb_by_dbcan"
DARKER_COLOR = (220/255, 220/255, 220/255)
LIGHT_COLOR = (240/255, 240/255, 240/255)



