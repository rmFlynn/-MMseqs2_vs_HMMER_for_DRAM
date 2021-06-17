"""Perform a HMMER scan"""
import os
import argparse
from mag_annotator.utils import run_process
from datetime import datetime

def working_dir_setup(working_dir):
    """Set up the dir"""
    assert not os.path.exists(working_dir), "Results already exist"
    os.mkdir(working_dir)

def main():
    """Maine function"""
    threads = 2
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genes_faa", type=str, default=None,
                        help="The protein data to be annotated usaly an fa")
    parser.add_argument("-m", "--hmm", type=str, default=None,
                        help="The HMM file")
    parser.add_argument("-o", "--output_dir", type=str, default=None,
                        help="The output directory")
    parser.add_argument("-n", "--name", type=str, default=None,
                        help="A name for the run")
    parser.add_argument("-c", "--clock_run", type=bool, default=False,
                        help="To time or not time the run, default not")
    par = parser.parse_args()

    working_dir = os.path.join(
        par.output_dir,
        par.name + "_results_"+ datetime.now().strftime("%Y_%m_%d_%H"))
    working_dir_setup(working_dir)
    os.chdir(working_dir)

    if par.clock_run :
        if not os.path.exists("run_times"):
            os.mkdir("run_times")
        time_string = [
            "/usr/bin/time", "-o",
            os.path.join(
                working_dir,
                "run_times",
                "time_hmmer"
            ),
            "-v",]
    else:
        time_string = []

    # now we can run
    run_process(time_string + [
        'hmmsearch',
        '--domtblout',
        os.path.join(
            working_dir,
            'hmmer_results.b6'
        ),
        '--cpu', str(threads),
        par.hmm,
        par.genes_faa],
        verbose=True)

if __name__ == "__main__":
    main()
