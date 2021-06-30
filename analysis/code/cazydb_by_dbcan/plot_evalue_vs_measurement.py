"""Download format data and make a plot"""
import sys
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from shared_tools import OUT_PLOT_PATH, OUT_DATA_PATH
sys.path.insert(0, "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/"
                "analysis/code/")
from output_tools.plot_tools import LIGHT_COLOR
from data_tools.mmseqs_vs_hmmer_table import split_out_eval_sens

data = pd.read_pickle(os.path.join(
    OUT_DATA_PATH, 'gold_standard_comp_stats.pkl'))
data.set_index('name', drop=True, inplace=True)
data = pd.concat([
    data[data.index.values == 'hmmer_default'],
    split_out_eval_sens(
        data[data.index.values != 'hmmer_default'])
])
# limit to 7.5 for mmseqs_sensitivity
data = data[(data['sensitivity'] == 7.5) | (data['sensitivity'].isna())]
data.reset_index(inplace=True)
data = data.melt(id_vars=['e-value', 'sensitivity', 'name'])
data = data[data['variable'] != 'True Negitive']
data = data.rename(columns={'variable': 'measure', 'value': 'count'})
data['e-value_exp'] = np.log10(data['e-value'])

pal = sns.cubehelix_palette(3,
                            light=0.7, start=.1,
                            rot=-1.75)
sns.set_style("whitegrid")
plt.clf()
fig, ax = plt.subplots(figsize=(6, 6))
sns.lineplot(style='measure',
             y="count",
             x='e-value_exp',
             hue='measure',
             palette=pal,
             ax=ax,
             data=data[data['name'] != 'hmmer_default'],
             markers=True
             )
ax.legend(title='MMseqs Sensitivity Setting')
for i, dat in data[data['name'] == 'hmmer_default'].reset_index().iterrows():
    ax.scatter(x=-15, y=dat['count'], color=pal[i])
    ax.annotate("HMMER3", xy=(-15, dat['count']),
                horizontalalignment='right')

ax.set_title("")
# axs[0].legend(loc='lower center', ncol=3, bbox_to_anchor=(3, -0.25))
ax.set_xlabel("MMseqs E-Value Exponent")

fig.subplots_adjust(wspace=0.01, bottom=0.2)
fig.patch.set_facecolor(LIGHT_COLOR)
fig.savefig(os.path.join(
    OUT_PLOT_PATH, 'mmseqs_vs_hmmer_gold_comp.png'))

def main():
    pass


if __name__ == '__main__':
    main()





