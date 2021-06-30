"""All functions and objects shared across projects anotated by dbCAN"""
import re
import numpy as np
import pandas as pd
from reading_output_tools import list_proteins

RAW_ALIGNMENTS = \
    "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/dbcan/dbCAN-fam-aln-V9.tar.gz"

PROTEIN_RE = r"^[A-Z]+\d+$"

def flatten_protein_ids(all_proteins):
    all_proteins = pd.Series(all_proteins).str.\
        split('|', expand = True)
    # NOTE: I checked the first column dose not match the regex
    all_proteins = all_proteins.loc[:, 1:].values.flatten()
    all_proteins = all_proteins[all_proteins != None]
    all_proteins = all_proteins[all_proteins != '']
    return all_proteins

def filter_protein_ids(proteins, protein_re):
    proteins = np.array([str(x).split('_')[0] for x in proteins])
    keep_proteins = np.array([x for x in  proteins
                              if re.match(protein_re, x)])
    drop_proteins = np.array([x for x in  proteins
                              if not re.match(protein_re, x)])
    return np.unique(keep_proteins), np.unique(drop_proteins)

def split_out_protein_list(proteins, protein_re):
    protein_df = pd.DataFrame(
        proteins.apply(
            lambda x: filter_protein_ids(x.split('|')[1:], protein_re)).to_list(),
        index=proteins.index, columns=['keep', 'drop'])
    drop = np.unique(np.concatenate(protein_df['drop'].values))
    return protein_df['keep'], drop

def make_mmseqs_vs_cazydb_truth(df, protein_re):
    data = df.copy()
    data['annotation'], _ = filter_protein_ids(data['annotation'], protein_re)
    data['true IDs'], _= split_out_protein_list(data['ProteinID'], protein_re)
    data['true'] = data.apply(
        lambda x: 1 if x['annotation'] in x['true IDs'] else 0, 1)
    return data

def count_proteins(protein_file, protein_re):
    return len(filter_protein_ids(flatten_protein_ids(list_proteins(protein_file)), protein_re)[0])

