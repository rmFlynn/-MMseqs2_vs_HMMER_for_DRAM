"""Make a plot showing the relationship of the hmmer and mmseqs e-values around the decision point"""
import os
import sys

from shared_tools import OUT_PLOT_PATH, MMSEQS_SWEEP_OUTPUT, MMSEQS_ARGS, \
    HMMER_ARGS
sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from output_tools.plot_evalue_relation import make_plot_data, make_scater


if __name__ == "__main__":
    data = make_plot_data(MMSEQS_SWEEP_OUTPUT, MMSEQS_ARGS,HMMER_ARGS)

    make_scater(data, OUT_PLOT_PATH, zoom=50)
    make_scater(data, OUT_PLOT_PATH, zoom=100)

