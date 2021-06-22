"""Download format data and make a plot"""
import sys
import os
import numpy as np
import pandas as pd
from shared_tools import parse_hmmsearch_domtblout, \
    read_mmseqs_results, list_vog_protines, list_vogs, \
    make_hmmer_vs_mmseqs_df, MMSEQS_SWEEP_OUTPUT
"""

this part of the analysis will use the default hmmer data and a representative
from the MMseqs sweep specifically the e-value = 15, sensitivity = 7.5 data
sets). These two were selected as the balanced option. It is suspected that
these settings will come the closest to matching, and so these settings will
be most representative of the final outcome. There is a risk however that
these more restrictive settings will also result in missing some anomalies
that are randomly filtered out of only these results. While the concern is
valid I have decided to take that risk so that the most important results will
be checked more throughly.

"""
hmmer_df = parse_hmmsearch_domtblout()
mmseqs_df = read_mmseqs_results(os.path.join(
    MMSEQS_SWEEP_OUTPUT, "evalue_1e-15_sens_7.500000"))
"""

Checking Protein IDs

The protein ides are different for each data set. Here we will do some simple
checks to see if they are the same as those in the protein file.

First we will check that all protein ids appear only once.
NOTE: This is specific to vogdb as we only take the top result here.

Demonstrare the validity of the aproche:

"""
invalid_tester = pd.DataFrame({
    'ProteinID': ['1000664.YP_009506390.1', '1000664.YP_009506390.1'],
})
valid_tester = pd.DataFrame({
    'ProteinID': ['1000664.YP_009506390.1', '1000664.YP_009506390.2'],
})
invalid_tester[invalid_tester['ProteinID'].duplicated()]
valid_tester[valid_tester['ProteinID'].duplicated()]
"""
Apply
"""
hmmer_df[hmmer_df['ProteinID'].duplicated()]
mmseqs_df[mmseqs_df['ProteinID'].duplicated()]
"""
There was no output, indicating no duplications.


Now we ask, are all the protein IDs in the original protein data base?
"""
proteins = pd.DataFrame({'ProteinID': list_vog_protines()})

def set_up_comparison(data, comp, on):
    data = data.copy()
    comp['is in in'] = 1
    data['is in out'] = 1
    comp_df = pd.merge(data, comp, on=on, how='outer')
    comp_df['is in out'].fillna(0, inplace=True)
    comp_df['is in in'].fillna(0, inplace=True)
    return comp_df

hmmer_protein = set_up_comparison(hmmer_df, proteins, on='ProteinID')
mmseqs_protein = set_up_comparison(mmseqs_df, proteins, on='ProteinID')
"""
Here we will list observations that are in the output but not in the input.
This would only happen if there is an error so we want no results.
"""
hmmer_protein[hmmer_protein['is in in'] == 0]
mmseqs_protein[mmseqs_protein['is in in'] == 0]
"it looks like we are OK on this for now. "

"""
Here we will list and count observations that are in the input but not in the
output. It is expected that some proteins will not get matched but we want it
to be few. and we don't want it to happen because of the wrong reasons.
"""
hmmer_protein[hmmer_protein['is in out'] == 0]
mmseqs_protein[mmseqs_protein['is in out'] == 0]
"""
there are 51480 observations for mmseqs and 49351 for hmmer. This is about
as we would expect.

Now we will check vogs in the same way.
"""
vogs = pd.DataFrame({'annotation': list_vogs()})
hmmer_vogs = set_up_comparison(hmmer_df, vogs, on='annotation')
mmseqs_vogs = set_up_comparison(mmseqs_df, vogs, on='annotation')
"""
Here we will list observations that are in the output but not in the input.
This would only happen if there is an error so we want no results.
"""
hmmer_vogs[hmmer_vogs['is in in'] == 0]
mmseqs_vogs[mmseqs_vogs['is in in'] == 0]
"""
it looks like we are OK on this for now.

Here we will list and count observations that are in the input but not in the
output. It is expected that some vogs will not get matched but we want it
to be few. and we don't want it to happen because of the wrong reasons.
"""
hmmer_vogs[hmmer_vogs['is in out'] == 0]
mmseqs_vogs[mmseqs_vogs['is in out'] == 0]
"""
there are 25683 observations for mmseqs and 25567 hmmer. This
is about as we would expect seeing as this is a small data set.


"""
mm_v_hm = make_hmmer_vs_mmseqs_df(hmmer_df, mmseqs_df)
mm_v_hm = mm_v_hm[(mm_v_hm['mmseqs_pred'] == 0) | \
                  (mm_v_hm['hmmer_pred'] == 0)]
mm_v_hm.sort_values("ProteinID")
mm_v_hm.sort_values("annotation")
# mm_v_hm.to_csv()









