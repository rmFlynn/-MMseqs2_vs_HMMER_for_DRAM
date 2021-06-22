"""This will make the stats table and tables for checking results"""
import sys
import os
from io import StringIO
from multiprocessing import Pool
from subprocess import Popen, PIPE
from sklearn.metrics import confusion_matrix
import pandas as pd
from shared_tools import parse_hmmsearch_domtblout, \
    read_mmseqs_results


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

def make_mmseqs_stats_df(mmseqs_result_path, vog_truth, MMSEQS_SWEEP_OUTPUT,
                         total_negitive):
    cfmx = mmseqs_cfmx(MMSEQS_SWEEP_OUTPUT + mmseqs_result_path, vog_truth,
                       total_negitive)
    return make_stats_df(mmseqs_result_path, cfmx)

def wrap_make_mmseqs_stats_df(args):
    return make_mmseqs_stats_df(*args)


def read_vog_truth(file_name:str):
    data = pd.read_csv(file_name, sep='\t')
    data = pd.concat(
        [data[data.columns[:-1]],
         data['ProteinIDs'].str.split(',', expand=True)],
        axis=1).\
        melt( id_vars=data.columns[:-1],
             var_name="drop",
             value_name='ProteinID').\
        drop(['drop'], axis=1).\
        dropna(subset=['ProteinID'])
    # data.drop_duplicates(subset=['ProteinID'], inplace=True)
    data.rename(columns={'#GroupName':'annotation'}, inplace=True)
    data = data[['ProteinID', 'annotation']]
    return data


def make_hmmmer_df(hmmer_results, vog_truth):
    hmmer = hmmer_results.copy()
    hmmer["hmmer_pred"] = 1
    vog_truth["true_pred"] = 1
    hmmer.rename(columns={"e-value":"hmmer_e-value"}, inplace=True)
    data = pd.merge(hmmer, vog_truth, how='outer',
                        on=['ProteinID', 'annotation'])
    data['hmmer_e-value'].fillna(10, inplace=True)
    data['hmmer_pred'].fillna(0, inplace=True)
    data['true_pred'].fillna(0, inplace=True)
    return data

def make_mmseqs_df(mmseqs_results, vog_truth):
    mmseqs = mmseqs_results.copy()
    mmseqs["mmseqs_pred"] = 1
    vog_truth["true_pred"] = 1
    mmseqs.rename(columns={"e-value":"mmseqs_e-value"}, inplace=True)
    data = pd.merge(mmseqs, vog_truth, how='outer',
                        on=['ProteinID', 'annotation'])
    data['mmseqs_e-value'].fillna(10, inplace=True)
    data['mmseqs_pred'].fillna(0, inplace=True)
    data['true_pred'].fillna(0, inplace=True)
    return data

def hmmer_cfmx(hmmer_results, vog_truth, total_negitive):
    base = make_hmmmer_df(hmmer_results, vog_truth)
    cfmx = confusion_matrix(base['true_pred'], base['hmmer_pred'])
    cfmx = make_cf_real(cfmx, total_negitive)
    return cfmx

def mmseqs_cfmx(mmseqs_path, vog_truth, total_negitive):
    mmseqs_results = read_mmseqs_results(mmseqs_path)
    mmseqs_cmp = make_mmseqs_df(mmseqs_results, vog_truth)
    cfmx = confusion_matrix(mmseqs_cmp['true_pred'],
                            mmseqs_cmp['mmseqs_pred'])
    cfmx = make_cf_real(cfmx, total_negitive)
    return cfmx

def make_final_df():
    """
    Here we are performing the gold standard analysis baseline analysis

    The Date is read in using functions from reading_output_tools in the root
    of the analysis folder. All folders in results should have only one file
    to avoid confusion, still the paths should be checked.
    """
    vog_truth = read_vog_truth()
    hmmer_results = parse_hmmsearch_domtblout()

    """
    The total negatives protein counts are made here. The process should be
    simple but unique to each gold standard analysis.
    """
    vog_count = count_vogs()
    protein_count = count_vog_proteins(RAW_PROTEINS
        )
    total_negitive = protein_count * vog_count - len(vog_truth)

    hmmer_stats = make_stats_df(
        'hmmer_default',
        hmmer_cfmx(hmmer_results, vog_truth, total_negitive))

    with Pool(10) as pool:
        mmseqs_dfs = pool.map(
            wrap_make_mmseqs_stats_df,
            [(i, vog_truth, MMSEQS_SWEEP_OUTPUT, total_negitive)
             for i in os.listdir(MMSEQS_SWEEP_OUTPUT)])

    mmseqs_dfs.append(hmmer_stats)
    all_stats = pd.concat(mmseqs_dfs)
    all_stats.sort_values('name', ascending=False, inplace=True)
    all_stats.reset_index(inplace=True, drop=True)
    all_stats.to_pickle(os.path.join([
        OUT_DATA_PATH, "all_confusion_matrix_and_stats.pkl"]))
    all_stats.to_csv(os.path.join([
        OUT_DATA_PATH, "all_confusion_matrix_and_stats.pkl"]), index=False)

if __name__ == '__main__':
    make_final_df()


