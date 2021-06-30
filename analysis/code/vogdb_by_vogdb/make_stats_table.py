"""This will make the stats table and tables for checking results"""
import sys
import os
from multiprocessing import Pool
from sklearn.metrics import confusion_matrix
import pandas as pd
from shared_tools import MMSEQS_SWEEP_OUTPUT, HMMER_ARGS, MMSEQS_ARGS, \
    VOG_TRUTH_FILE, RAW_PROTEINS, OUT_DATA_PATH
sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from reading_output_tools import parse_hmmsearch_domtblout, \
    read_mmseqs_results, set_up_comparison, list_proteins, \
    list_annotations, make_hmmer_vs_mmseqs_df
from vog_tools import RAW_ALIGNMENTS, read_vog_truth, \
    make_hmmmer_vs_vog_truth_df, make_mmseqs_vs_vog_truth_df

"""Tools to make the stats using modified confusion_matrixis"""

def make_cf_real(cfmx, total_negitive):
    tn, fp, fn, tp = cfmx.ravel() # safe way to do this
    cfmx[0,0] = total_negitive - fp
    return cfmx


def sensitivity(cfmx):
    tn, fp, fn, tp = cfmx.ravel() # safe way to do this
    return tp / (tp + fn)

def specificity(cfmx):
    tn, fp, fn, tp = cfmx.ravel() # safe way to do this
    return tn / (fp + tn)

def full_stats(cfmx):
    sp = specificity(cfmx)
    se = sensitivity(cfmx)
    ba = (sp + se) / 2
    return sp, se, ba

def make_stats_df(name, cfmx):
    sp, se, ba = full_stats(cfmx)
    df = pd.DataFrame({
        'name': name,
        'confusion_matrix': [cfmx],
        'specificity': sp,
        'sensitivity': se,
        'balanced_accuracy': ba
    })
    return df

"""end"""

def make_mmseqs_stats_df(mmseqs_result_path, vog_truth, mmseqs_sweep_output,
                         total_negitive, mmseqs_args):
    cfmx = mmseqs_cfmx(mmseqs_sweep_output + mmseqs_result_path, vog_truth,
                       total_negitive, mmseqs_args)
    return make_stats_df(mmseqs_result_path, cfmx)

def wrap_make_mmseqs_stats_df(args):
    return make_mmseqs_stats_df(*args)

def count_vog_proteins(raw_proteins):
    """
    :param file_name: The name of the protein fa file
    :returns: a count of all proteins
    """
    return len(list_proteins(raw_proteins))

def count_vogs(raw_alignments):
    """
    :param file_name: The raw aliments file
    :returns: the count of all vogs
    """
    return len(list_annotations(raw_alignments, ['.msa']))

def hmmer_cfmx(hmmer_results, vog_truth, total_negitive):
    base = make_hmmmer_vs_vog_truth_df(hmmer_results, vog_truth)
    cfmx = confusion_matrix(base['true_pred'], base['hmmer_pred'])
    cfmx = make_cf_real(cfmx, total_negitive)
    return cfmx

def mmseqs_cfmx(mmseqs_path, vog_truth, total_negitive, mmseqs_args):
    mmseqs_results = read_mmseqs_results(mmseqs_path, **mmseqs_args)
    mmseqs_cmp = make_mmseqs_vs_vog_truth_df(mmseqs_results, vog_truth)
    cfmx = confusion_matrix(mmseqs_cmp['true_pred'],
                            mmseqs_cmp['mmseqs_pred'])
    cfmx = make_cf_real(cfmx, total_negitive)
    return cfmx

def make_final_df(hmmer_args, mmseqs_args, mmseqs_sweep_output, raw_proteins,
                  raw_alignments, vog_truth_file, out_data_path):
    """
    Here we are performing the gold standard analysis baseline analysis

    The Date is read in using functions from reading_output_tools in the root
    of the analysis folder. All folders in results should have only one file
    to avoid confusion, still the paths should be checked.
    """
    vog_truth = read_vog_truth(vog_truth_file)
    hmmer_results = parse_hmmsearch_domtblout(**hmmer_args)
    """
    The total negatives protein counts are made here. The process should be
    simple but unique to each gold standard analysis.
    """
    alignment_count = count_vogs(raw_alignments)
    protein_count = count_vog_proteins(raw_proteins)
    total_negitive = protein_count * alignment_count - len(vog_truth)

    hmmer_stats = make_stats_df(
        'hmmer_default',
        hmmer_cfmx(hmmer_results, vog_truth, total_negitive))

    with Pool(10) as pool:
        mmseqs_dfs = pool.map(
            wrap_make_mmseqs_stats_df,
            [(i, vog_truth, mmseqs_sweep_output, total_negitive, mmseqs_args)
             for i in os.listdir(MMSEQS_SWEEP_OUTPUT)])

    mmseqs_dfs.append(hmmer_stats)
    all_stats = pd.concat(mmseqs_dfs)
    all_stats.sort_values('name', ascending=False, inplace=True)
    all_stats.reset_index(inplace=True, drop=True)
    all_stats.to_csv(os.path.join(
        out_data_path, 'gold_standard_comp_stats.csv'), index=False)
    all_stats.to_pickle(os.path.join(
        out_data_path, 'gold_standard_comp_stats.pkl'))

if __name__ == '__main__':
    make_final_df(HMMER_ARGS, MMSEQS_ARGS, MMSEQS_SWEEP_OUTPUT, RAW_PROTEINS,
                  RAW_ALIGNMENTS, VOG_TRUTH_FILE, OUT_DATA_PATH)

