"""Make latex tabels that show diferance"""
import pandas as pd
import numpy as np

def color_cfmx(cfmx):
    cfmxs = np.chararray((2, 2), 50)
    cfmxs[0,1] = '{\color{red}' + str(cfmx[0,1]) + '}' \
        if cfmx[0,1] > 0 else '{\color{green}' + str(cfmx[0,1]) + '}'
    cfmxs[1,0] = '{\color{red}' + str(cfmx[1,0]) + '}' \
        if cfmx[1,0] > 0 else '{\color{green}' + str(cfmx[1,0]) + '}'
    cfmxs[0,0] = '{\color{green}' + str(cfmx[0,0]) + '}' \
        if cfmx[0,0] > 0 else '{\color{red}' + str(cfmx[0,0]) + '}'
    cfmxs[1,1] = '{\color{green}' + str(cfmx[1,1]) + '}' \
        if cfmx[1,1] > 0 else '{\color{red}' + str(cfmx[1,1]) + '}'
    return cfmxs

def main():
    all_stats = pd.concat([
        pd.read_pickle("all_stats.pkl"),
        pd.read_pickle("timed_stats.pkl")])
    all_stats.set_index('name', inplace=True)
    dif_df = all_stats.loc[[
        'evalue_1e-14_sens_7.500000',
        'evalue_1e-15_sens_7.500000',
        'evalue_1e-16_sens_7.500000', ]]
    hmmer = all_stats.loc['hmmer_default']['confusion_matrix'].values[0]
    cfs = dif_df['confusion_matrix'].apply(lambda x: x - hmmer)

    for i, j in pd.DataFrame(cfs).iterrows():
        with open(os.path.join(
            OUT_LATEX_PATH, "dif_" + i +".tex"), "w") as f:
            f.write(
                pd.DataFrame(
                    color_cfmx(j[0]),
                    columns=['0', '1'],
                    index=['0', '1']
                ).to_latex(escape=False)
            )

if __name__ == '__main__':
    main()
