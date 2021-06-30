"""This will make the stats table and tables for checking results"""
import sys
import os
from multiprocessing import Pool
from sklearn.metrics import confusion_matrix
import pandas as pd
import numpy as np
from shared_tools import MMSEQS_SWEEP_OUTPUT, HMMER_ARGS, MMSEQS_ARGS, \
    RAW_PROTEINS, OUT_DATA_PATH, PROTEIN_RE
sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from reading_output_tools import parse_hmmsearch_domtblout, \
    read_mmseqs_results, set_up_comparison, list_proteins, \
    list_annotations, make_hmmer_vs_mmseqs_df
from dbcan_tools import RAW_ALIGNMENTS, split_out_protein_list

def sensitivity(cfmx):
    tn, fp, fn, tp = cfmx.ravel() # safe way to do this
    return tp / (tp + fn)

def specificity(cfmx):
    tn, fp, fn, tp = cfmx.ravel() # safe way to do this
    return tn / (fp + tn)

def make_stats_line(name, cfmx):
    tn, fp, fn, tp = cfmx.ravel() # safe way to do this
    df = pd.DataFrame({
        'name': name,
        'False Positve': fp,
        'False Negitive': fn,
        'True Positive': tp,
        'True Negitive': tn
    }, index=name)
    return df

def make_mmseqs_stats_line(mmseqs_result_path, mmseqs_sweep_output,
                           mmseqs_args, gold_truth):
    results = read_mmseqs_results(os.path.join(
        mmseqs_sweep_output, mmseqs_result_path), **mmseqs_args)
    cfmx = make_cfmx(results, gold_truth)
    return make_stats_line(mmseqs_result_path, cfmx)
def wrap_make_mmseqs_stats_line(args):
    return make_mmseqs_stats_line(*args)

def make_cfmx(in_results, in_gold_truth):
    results = in_results.copy()
    gold_truth = in_gold_truth.copy()
    gold_truth['true'] = 1
    results['pred'] = 1
    stats = pd.merge(results, gold_truth,
             on=['ProteinID', 'annotation'],  how='outer')
    stats[['pred', 'true']] = stats[['pred', 'true']].fillna(0)
    cfmx = confusion_matrix(stats['true'], stats['pred'])
    return cfmx

# results['true_ids'], _ = split_out_protein_list(results['ProteinID'],
#                                                 protein_re)
# results['true_pred'] = results.apply(
#     lambda x: int(x['annotation'] in x['true_ids']), 1)
# cfmx = confusion_matrix(results['true_pred'], np.ones(len(results)))

def read_make_gold_truth(raw_proteins, protein_re):
    gold_truth = pd.DataFrame({'ProteinID': list_proteins(RAW_PROTEINS)})
    gold_truth = pd.concat([
        gold_truth[['ProteinID']],
        split_out_protein_list(
            gold_truth['ProteinID'],
            protein_re
        )[0].apply('|'.join).str.split('|', expand=True)
    ], 1).melt(id_vars='ProteinID').dropna()
    gold_truth = gold_truth[['ProteinID', 'value']].drop_duplicates()
    # gold_truth[gold_truth['ProteinID'] == \
    #            'QAT16263.1|GH133|GH13|GH23|GH2|GT2|GT4|GT5|']
    gold_truth.columns = ['ProteinID', 'annotation']
    return gold_truth


def make_final_df(hmmer_args, mmseqs_args, mmseqs_sweep_output, raw_proteins,
                  raw_alignments, out_data_path, protein_re):
    """
    Here we are performing the gold standard analysis baseline analysis

    """
    gold_truth = read_make_gold_truth(raw_proteins, protein_re)
    hmmer_results = parse_hmmsearch_domtblout(**hmmer_args)


    hmmer_stats = make_stats_line(
        'hmmer_default',
        make_cfmx(hmmer_results, gold_truth))

    with Pool(10) as pool:
        mmseqs_dfs = pool.map(
            wrap_make_mmseqs_stats_line,
            [(i, mmseqs_sweep_output, mmseqs_args, gold_truth)
             for i in os.listdir(mmseqs_sweep_output)])

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
                  RAW_ALIGNMENTS, OUT_DATA_PATH, PROTEIN_RE)

