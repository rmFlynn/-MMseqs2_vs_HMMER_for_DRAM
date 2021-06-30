"""Make tables that compare the results of mmseqs and vogdb"""
import sys
from shared_tools import MMSEQS_SWEEP_OUTPUT, HMMER_ARGS, MMSEQS_ARGS, \
    OUT_DATA_PATH
sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from data_tools.mmseqs_vs_hmmer_table import make_count_table

if __name__ == '__main__':
    make_count_table(MMSEQS_SWEEP_OUTPUT, OUT_DATA_PATH, HMMER_ARGS, MMSEQS_ARGS)

# import os
# os.system('cat ../code_tests/hmmer_test_file')
# HMMER_ARGS['file_name'] = '../code_tests/hmmer_test_file'
# os.system('mkdir ../code_tests')
# out = parse_hmmsearch_domtblout(**HMMER_ARGS)
# results = out
# results = results[results['annotation'].apply( lambda x: bool(re.match(HMMER_ARGS['filter_regex'], x)))]


