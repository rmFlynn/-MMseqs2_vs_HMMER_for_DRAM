"""
This is the first file that should be run and requiers actual thoght
"""
import sys
import numpy as np
import pandas as pd
from shared_tools import parse_hmmsearch_domtblout, \
    read_mmseqs_results, filter_protein_ids, split_out_protein_list, \
    MMSEQS_SWEEP_OUTPUT, flatten_protein_ids, list_anotations

hmmer_df = parse_hmmsearch_domtblout()

mmseqs_df = read_mmseqs_results(os.path.join(
    MMSEQS_SWEEP_OUTPUT, "evalue_1e-15_sens_7.500000"))

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

The first step is specific to cazydb an dbcan. We filter the protein ids with
regex comand wich is hard to trust. Fortunatly it is easy to check just look
the output of the falowing comands. Good should have valid protein ids and bad should not.

The best place to start is with the faa file itself.
"""
good, bad = filter_protein_ids(flatten_protein_ids(list_proteins()))
good
bad
"""
This seam to have worked well
"""
good, bad = filter_protein_ids(mmseqs_df['annotation'].values)
good
bad
"""
It may take some time. Note that 'CBM35inCE17', 'dockerin' are not bad
proteins but are not in the target data base.
"""
good, bad = filter_protein_ids(hmmer_df['annotation'].values)
good
bad
"""
Same story with 'SLH' and 'cohesin'
"""
good, bad = split_out_protein_list(mmseqs_df['ProteinID'].values)
np.unique(np.concatenate(good.values))
bad
"""
Again no complaints
"""
good, bad = split_out_protein_list(hmmer_df['ProteinID'].values)
np.unique(np.concatenate(good.values))
bad
"""
seems we are good here.

## Checking Protein IDs

The protein ides are different for each data set. Here we will do some simper
checks to see if they are the same as those in the protein file.

First we will check that all protein ids appear only once.
NOTE: This is specific to vogdb as we only take the top result here.

Demonstrare the validity of the aproche:

invalid_tester = pd.DataFrame({
    'ProteinID': ['1000664.YP_009506390.1', '1000664.YP_009506390.1'],
})
valid_tester = pd.DataFrame({
    'ProteinID': ['1000664.YP_009506390.1', '1000664.YP_009506390.2'],
})
invalid_tester[invalid_tester['ProteinID'].duplicated()]
valid_tester[valid_tester['ProteinID'].duplicated()]

Apply

hmmer_df[hmmer_df['ProteinID'].duplicated()]
mmseqs_df[mmseqs_df['ProteinID'].duplicated()]



Now we ask, are all the protein IDs in the original protein data base?
"""
proteins = pd.DataFrame({'ProteinID': list_proteins()})

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
"""
NOTE "Q9Y2B1.1" is interesting

Here we will list and count observations that are in the input but not in the
output. It is expected that some proteins will not get matched but we want it
to be few. and we don't want it to happen because of the wrong reasons.
"""
hmmer_protein[hmmer_protein['is in out'] == 0]
mmseqs_protein[mmseqs_protein['is in out'] == 0]
274559 / 1716043
639476 / 1716043
"""
there are 639476 (38%) observations for mmseqs and for 274559 (15%) hmmer. This is more
more than expected for mmseqs but it may be realated to the odd ides above

Now we will check annotations in the same whay.
"""
annotations = pd.DataFrame({
    'annotation':
    filter_protein_ids(list_anotations(['dbCAN-fam-aln/', '.aln']))[0]
})
hmmer_df['annotation'] = \
    hmmer_df['annotation'].str.split('_', expand=True)[0]
mmseqs_df['annotation'] = \
    mmseqs_df['annotation'].str.split('_', expand=True)[0]
hmmer_ano = set_up_comparison(hmmer_df, annotations, on='annotation')
mmseqs_ano = set_up_comparison(mmseqs_df, annotations, on='annotation')
"""
Here we will list observations that are in the output but not in the input.
This would hapen in dbCAN unlike in vog db, because these observations can't
be confermed or denied they are filterd out of the final analysis. Still we
don't want to many of these.
"""
hmmer_ano[hmmer_ano['is in in'] == 0]
mmseqs_ano[mmseqs_ano['is in in'] == 0]
"""
hmmer has many more out of scope annotation compared to mmseqs. It is also
clear that CBM35inCE17 sould posibly be included that can be fixed later.

NOTE: I cant garenty that removing these observations will not change the
result if there is a corlation with the models and the corectness.

Here we will list and count observations that are in the input but not in the
output. It is expected that some annotations will not get matched but we want
it to be few. and we don't want it to happen because of the wrong reasons.
"""
hmmer_ano[hmmer_ano['is in out'] == 0]
mmseqs_ano[mmseqs_ano['is in out'] == 0]

"""
Finally it is good to take a look at the default file reads
"""





