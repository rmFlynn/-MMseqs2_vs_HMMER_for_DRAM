"""Download format data and make a plot"""
import os
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE
from io import StringIO
from functools import reduce
from matplotlib import pyplot as plt
from sklearn.preprocessing import label_binarize
from sklearn import svm, datasets
from sklearn.metrics import confusion_matrix
from time import sleep
from sklearn.preprocessing import MinMaxScaler
from decimal import Decimal
from multiprocessing import Pool

BOUTFMT6_COLUMNS_MMSEQS = ['ProteinID', 'vog', 'seqIdentity', 'alnLen',
                        'mismatchCnt', 'gapOpenCnt', 'qStart', 'qEnd',
                        'tStart', 'tEnd', 'e-value', 'bitScore']
HMMSCAN_ALL_COLUMNS = ['query_id', 'query_ascession', 'query_length', 'target_id',
                       'target_ascession', 'target_length', 'full_evalue', 'full_score', 'full_bias', 'domain_number', 'domain_count', 'domain_cevalue', 'domain_ievalue', 'domain_score', 'domain_bias', 'target_start', 'target_end', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'accuracy', 'description']
HMMSCAN_COLUMN_TYPES = [str, str, int, str, str, int, float, float, float, int, int, float, float, float, float, int, int, int, int, int, int, float, str]

def read_vog_truth(path):
    data = pd.read_csv(path, sep='\t')
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
    data.rename(columns={'#GroupName':'vog'}, inplace=True)
    data = data[["ProteinID", "vog"]]
    return data

def parse_hmmsearch_domtblout(file_name:str, pct_cover=0.35):
    df_lines = list()
    for line in open(file_name):
        if not line.startswith('#'):
            line = line.split()
            line = line[:22] + [' '.join(line[22:])]
            df_lines.append(line)
    results = pd.DataFrame(df_lines, columns=HMMSCAN_ALL_COLUMNS)
    results[["target_start", "target_end", "target_length"]] = \
        results[["target_start", "target_end", "target_length"]].astype(float)
    results['pct_cover'] = results.apply(
        lambda x: (x['target_end'] - x['target_start']) / x['target_length'],
        axis=1)
    results = results[results['pct_cover'] > pct_cover]
    results = results[['query_id', 'target_id', 'full_evalue']]
    results.columns = ['ProteinID', 'vog', 'e-value']
    results = results.astype({'ProteinID': str, 'vog': str, 'e-value':float})
    # filter evalue to match curent
    results = results[results['e-value'] < 1e-15]
    return results

def count_proteins():
    all_proteins = pd.read_csv(
        StringIO(
            Popen(
                ["cat", "../../data/vog.proteins.all.20210525.fa",],
                stdout=PIPE).communicate()[0].decode('utf-8')
        ),
        header=None,
    )
    all_proteins.columns = ['proteins']
    all_proteins = all_proteins[
        all_proteins['proteins']
        .str.contains(">", na=False)]
    all_proteins['proteins'] = all_proteins['proteins'].str.replace(">", "")
    return len(all_proteins['proteins'].unique())

def count_vogs():
    all_vogs = pd.read_csv(
        StringIO(
            Popen(
                ["tar", "-tvf", "../../data/vog.raw_algs.20210525.tar.gz",],
                stdout=PIPE).communicate()[0].decode('utf-8')
        ),
        header=None,
    delimiter=r"\s+"
    )
    all_vogs =  all_vogs[[5]]
    all_vogs.columns = ['vog']
    all_vogs['vog'] = all_vogs['vog'].str.replace(".msa", "")
    return len(all_vogs['vog'].unique())

def make_hmmmer_df(hmmer_results, vog_truth):
    hmmer = hmmer_results.copy()
    hmmer["hmmer_pred"] = 1
    vog_truth["true_pred"] = 1
    hmmer.rename(columns={"e-value":"hmmer_e-value"}, inplace=True)
    data = pd.merge(hmmer, vog_truth, how='outer',
                        on=['ProteinID', 'vog'])
    data['hmmer_e-value'].fillna(10, inplace=True)
    data['hmmer_pred'].fillna(0, inplace=True)
    data['true_pred'].fillna(0, inplace=True)
    return data

def read_mmseqs_results(path, pct_cover=0.35):
    results = pd.read_csv(path, sep='\t',
                          header=None, names=BOUTFMT6_COLUMNS_MMSEQS)
    results['pct_cover'] = results.apply(
        lambda x: (x['tEnd'] - x['tStart']) / x['alnLen'],
        axis=1)
    results = results[results['pct_cover'] > pct_cover]
    results['vog'] = results['vog'].\
        str.replace('.msa', '', regex=False)
    results = results[['ProteinID', 'vog', 'e-value']]
    return results

def make_mmseqs_df(mmseqs_results, vog_truth):
    mmseqs = mmseqs_results.copy()
    mmseqs["mmseqs_pred"] = 1
    vog_truth["true_pred"] = 1
    mmseqs.rename(columns={"e-value":"mmseqs_e-value"}, inplace=True)
    data = pd.merge(mmseqs, vog_truth, how='outer',
                        on=['ProteinID', 'vog'])
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

def make_mmseqs_stats_df(mmseqs_result_path, vog_truth, sweep_output,
                         total_negitive):
    cfmx = mmseqs_cfmx(sweep_output + mmseqs_result_path, vog_truth,
                       total_negitive)
    return make_stats_df(mmseqs_result_path, cfmx)

def make_latex_stat_rep(df, f):
    f.write(df[['name', 'specificity', 'sensitivity', 'balanced_accuracy']].\
            to_latex(index_names=False, header=False))

def make_latex_cfmx_rep(df, f):
    # f.write('confusion_matrix:\n\n')
    f.write(pd.DataFrame(df['confusion_matrix'],
            columns=[0,1], index=[0,1]).to_latex(index=False, index_names=False, header=False))

def wrap_make_mmseqs_stats_df(args):
    return make_mmseqs_stats_df(*args)

def make_final_df():
    vog_truth = read_vog_truth("../../data/vog.members.20210525.tsv.gz")
    hmmer_results = parse_hmmsearch_domtblout("../../hmmer/hmmer_results.b6")
    sweep_output = "../../mmseqs/results_2021_06_03_11/sweep_output/"
    total_negitive = count_proteins() * count_vogs() - len(vog_truth)

    hmmer_stats = make_stats_df(
        'hmmer_default',
        hmmer_cfmx(hmmer_results, vog_truth, total_negitive))

    with Pool(10) as pool:
        mmseqs_dfs = pool.map(
            wrap_make_mmseqs_stats_df,
            [(i, vog_truth, sweep_output, total_negitive)
             for i in os.listdir(sweep_output)])

    mmseqs_dfs.append(hmmer_stats)
    all_stats = pd.concat(mmseqs_dfs)
    all_stats.sort_values('name', ascending=False, inplace=True)
    all_stats.reset_index(inplace=True, drop=True)
    all_stats.to_pickle("timed_stats.pkl")
    all_stats.to_csv("timed_stats.csv", index=False)

def make_select_tex():
    all_stats = pd.read_pickle("all_stats.pkl")
    with open("latex/best_balanced_accuracy_stats.tex", "w") as f:
        make_latex_stat_rep(
            all_stats.sort_values('balanced_accuracy').iloc[-1], f)

    with open("latex/best_balanced_accuracy_cfmx.tex", "w") as f:
        make_latex_cfmx_rep(
            all_stats.sort_values('balanced_accuracy').iloc[-1], f)

    with open("latex/best_sensitivity_stats.tex", "w") as f:
        make_latex_stat_rep(
            all_stats.sort_values('sensitivity').iloc[-1], f)

    with open("latex/best_sensitivity_cfmx.tex", "w") as f:
        make_latex_cfmx_rep(
            all_stats.sort_values('sensitivity').iloc[-1], f)

    with open("latex/best_specificity_stats.tex", "w") as f:
        make_latex_stat_rep(
            all_stats.sort_values('specificity').iloc[-1], f)

    with open("latex/best_specificity_cfmx.tex", "w") as f:
        make_latex_cfmx_rep(
            all_stats.sort_values('specificity').iloc[-1], f)

    with open("latex/hmmer_default_stats.tex", "w") as f:
        make_latex_stat_rep(
            all_stats[all_stats['name'] == 'hmmer_default'].iloc[0], f)

    with open("latex/hmmer_default_cfmx.tex", "w") as f:
        make_latex_cfmx_rep(
            all_stats[all_stats['name'] == 'hmmer_default'].iloc[0], f)

    all_stats['fp'] = all_stats['confusion_matrix'].apply(lambda x: x[0,1])
    all_stats['tp'] = all_stats['confusion_matrix'].apply(lambda x: x[1,1])
    hmmer_cfmx = all_stats[all_stats['name'] == 'hmmer_default']\
        ['confusion_matrix'].values[0]
    all_stats['dif'] = all_stats['confusion_matrix'].\
        apply(lambda x: np.sum(np.abs(x - hmmer_cfmx)))

    all_stats['tp'] = all_stats['confusion_matrix'].apply(lambda x: x[1,1])

    most_similar = all_stats.sort_values('dif').iloc[1]

    best_tp_better_fp = all_stats[all_stats['fp'] < \
              all_stats[all_stats['name'] == 'hmmer_default']['fp'].values[0]]\
        .sort_values('tp').iloc[-1]

    with open("latex/most_similar_stats.tex", "w") as f:
        make_latex_stat_rep(most_similar, f)

    with open("latex/most_similar_cfmx.tex", "w") as f:
        make_latex_cfmx_rep(most_similar, f)

    with open("latex/best_tp_less_fp_stats.tex", "w") as f:
        make_latex_stat_rep(best_tp_better_fp, f)

    with open("latex/best_tp_less_fp_cfmx.tex", "w") as f:
        make_latex_cfmx_rep(best_tp_better_fp, f)



def make_all_tex():
    all_stats = pd.read_pickle("timed_stats.pkl")
    for i in all_stats['name'].values:
        with open("latex/time_" + i + "_stats.tex", "w") as f:
            make_latex_stat_rep(
                all_stats[all_stats['name'] == i].iloc[0], f)

        with open("latex/time_" + i + "_cfmx.tex", "w") as f:
            make_latex_cfmx_rep(
                all_stats[all_stats['name'] == i].iloc[0], f)

if __name__ == '__main__':
    make_final_df()
    make_all_tex()





