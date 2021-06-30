"""Test the reading tools"""
import sys
import pytest
import pandas as pd
sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from reading_output_tools import parse_hmmsearch_domtblout, \
    read_mmseqs_results
from dbcan_tools import PROTEIN_RE

@pytest.mark.read_parse
def hmmer_un_arg():
    pass

import os
# os.system("cp .//hmmer_test_file .//hmmer_strip_arg_test_file")

@pytest.mark.read_parse
def hmmer_strip_arg():
    data1 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_strip_arg_test_file',
            strip_from_annotations=[], take_only_one=False)
    data2 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_strip_arg_test_file',
            strip_from_annotations=['.hmm', '.thistoo'], take_only_one=False)
    assert list(data1['annotation'].values) == \
           ['CBM10.hmm', 'CBM10.thistoo', 'CBM10hmm'], \
           'strip_from_annotations error'
    assert list(data2['annotation'].values) == \
                ['CBM10', 'CBM10hmm'], \
           'strip_from_annotations error'

# os.system("cp ./test_data/hmmer_test_file "
#           "./test_data/hmmer_take_one_test_file")
@pytest.mark.read_parse
def hmmer_take_one_arg():
    data1 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_take_one_test_file',
            strip_from_annotations=[], take_only_one=False)
    data2 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_take_one_test_file',
            strip_from_annotations=[], take_only_one=True)
    assert list(data1['annotation'].values) == ['justThis', 'notThis'], \
        'take_only_one error'
    assert list(data2['annotation'].values) == ['justThis'], \
        'take_only_one error'


@pytest.mark.read_parse
def hmmer_dbcan_under_arg():
    data1 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_dbcan_under_test_file',
            strip_from_annotations=[], take_only_one=False)
    data2 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_dbcan_under_test_file',
            strip_from_annotations=[], take_only_one=True,
            filter_dbcan_under=True)
    assert len(data1) == 2, 'filter_dbcan_under error'
    assert list(data2['annotation'].values) == ['CBM10'], \
        'filter_dbcan_under error'

@pytest.mark.read_parse
def hmmer_regex_arg():
    data1 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_regex_test_file',
            strip_from_annotations=['.hmm'], take_only_one=False,
            filter_dbcan_under=True)
    data2 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_regex_test_file',
            strip_from_annotations=['.hmm'], take_only_one=False,
            filter_dbcan_under=True, filter_regex=PROTEIN_RE)
    assert len(data1) == 7, 'take_only_one error'
    assert list(data2['annotation'].values) == ['CBM10'], \
        'remove_dup_pairs error'

@pytest.mark.read_parse
def hmmer_limit_e_arg():
    data1 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_limit_e_test_file',
            strip_from_annotations=['.hmm'], take_only_one=False, limit_e=1e-2)
    data2 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_limit_e_test_file',
            strip_from_annotations=['.hmm'], take_only_one=False)
    assert len(data1) == 7, 'take_only_one error'
    assert len(data2) == 1, 'take_only_one error'

@pytest.mark.read_parse
def hmmer_pct_cover_arg():
    pass

@pytest.mark.read_parse
def hmmer_remove_dup_arg():
    data1 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_remove_dup_test_file',
            strip_from_annotations=[], take_only_one=False,
            remove_dup_pairs=False)
    data2 = parse_hmmsearch_domtblout(
            file_name='./test_data/hmmer_remove_dup_test_file',
            strip_from_annotations=[], take_only_one=False)
    assert len(data1) == 6, 'remove_dup_pairs take_only_one error'
    assert list(data2['annotation'].values) == ['only1', 'only2', 'only2'], \
        'remove_dup_pairs error'

@pytest.mark.read_parse
def mmseqs_remove_dup_arg():
    pass

@pytest.mark.read_parse
def mmseqs_un_arg():
    pass

