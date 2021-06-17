#!/home/rmflynn/miniconda3/envs/DRAM/bin/python
"""Download, format data and make a sweep of mmseqs prams"""

import os
import time
from subprocess import Popen
from datetime import datetime
import argparse

BOUTFMT6_COLUMNS = ['qId', 'tId', 'seqIdentity', 'alnLen', 'mismatchCnt',
                    'gapOpenCnt', 'qStart', 'qEnd', 'tStart', 'tEnd', 'eVal',
                    'bitScore']

###############
# Vogdb setup #
###############

def working_dir_setup(working_dir):
    """Set up the dir"""
    assert not os.path.exists(working_dir), "Results already exist"
    os.mkdir(working_dir)

def process_profile(working_dir, threads:int, raw_aln:str):
    os.chdir(working_dir)
    os.mkdir("mmseqs_profile")
    os.chdir("mmseqs_profile")
    Popen(['mmseqs', 'tar2db',
           raw_aln, # Input
           "vog_msa",
           # Output
           "--output-dbtype", "11",
           # "--tar-include", "'.*msa$'",
           "--threads",  str(threads)]).wait()
    os.system("mmseqs msa2profile vog_msa profile --match-mode 1")
    os.chdir(working_dir)

def process_target(working_dir, threads:int, protein_data:str, sensitivity:float=None):
    """
    process proteins to be annotated

    :param working_dir:
    :param threads:
    :param protein_data:
    :param sensitivity:
    """
    os.chdir(working_dir)
    os.system("mkdir to_annotate")
    Popen(['mmseqs', 'createdb',
          protein_data,
          'to_annotate/vog.proteins.all']).wait()
    # indexes DB
    if sensitivity is not None:
        os.system("mkdir to_annotate/temp")
        Popen(['mmseqs', 'createindex', 'to_annotate/vog.proteins.all',
               'to_annotate/temp', '-k', '6', '-s', str(sensitivity),
               '--threads', str(threads)]).wait()

def mmseqs_search(working_dir, sensitivity, evalue, threads, clock_run):
    """
    search against the db

    :param working_dir:
    :param sensitivity:
    :param evalue:
    :param threads:
    :param clock_run:
    """
    os.chdir(working_dir)
    if not os.path.exists("raw_sweep_output"):
        os.mkdir("raw_sweep_output")
    if not os.path.exists("sweep_output"):
        os.mkdir("sweep_output")
    if clock_run :
        if not os.path.exists("run_times"):
            os.mkdir("run_times")
        time_string = [
            "/usr/bin/time", "-o",
            os.path.join(
                working_dir,
                "run_times",
                "time_mmseqs_{:e}_{:.1f}_{}".format(evalue, sensitivity, threads)
            ),
            "-v",]
    else:
        time_string = []
    Popen(time_string + ['mmseqs', 'search',
           'to_annotate/vog.proteins.all',
           'mmseqs_profile/profile',
           'raw_sweep_output/evalue_%s_sens_%f' %(str(evalue), sensitivity),
           'temp', '--threads', str(threads), '-e', str(evalue), '-s',
           str(sensitivity), '-k', '6']).wait()
    Popen(['mmseqs', 'convertalis',
           'to_annotate/vog.proteins.all',
           'mmseqs_profile/profile',
           'raw_sweep_output/evalue_%s_sens_%f' %(str(evalue), sensitivity),
           'sweep_output/evalue_%s_sens_%f' %(str(evalue), sensitivity),
           ]).wait()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--threads", type=int, default=2,
                        help="number of threads for mmseqs")
    parser.add_argument("-p", "--protein_data", type=str, default=None,
                        help="The protein data to be annotated usaly an fa")
    parser.add_argument("-r", "--raw_aln", type=str, default=None,
                        help="The raw sequence alinemnts in a tar file")
    parser.add_argument("-n", "--name", type=str, default=None,
                        help="A name for the run")
    parser.add_argument("-o", "--output_dir", type=str, default=None,
                        help="The output directory")
    parser.add_argument("-c", "--clock_run", type=bool, default=False,
                        help="To time or not time the run, default not")
    par = parser.parse_args()
    working_dir = os.path.join(
        output_dir,
        par.name + "_results_"+ datetime.now().strftime("%Y_%m_%d_%H"))
    par = parser.parse_args()
    working_dir_setup(working_dir)
    process_profile(working_dir, par.threads, par.raw_aln)
    process_protine(working_dir, par.threads, par.protein_data)

    for sens in range(12, 16, 1):
        for evalue in list(range(-20, 0, 2)) + [-15]:
            mmseqs_search(working_dir, sens / 2, float("1e%i"%evalue), par.threads, par.clock_run)


