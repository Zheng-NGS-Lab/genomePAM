import glob
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import skew
import re
import numpy as np

def generate_dataframe(IdentifiedOfftargetsFile,runDate):
    LibID=IdentifiedOfftargetsFile.split('/')[-1].split('_')[0]
    d1 = pd.read_csv(IdentifiedOfftargetsFile, sep='\t', header=0, dtype=str, na_filter=False)
    target = d1['TargetSequence'][0]
    spacer = re.sub('N{2,}', '',target)
    PAMcand = target.replace(spacer, '')
    lenPAM = len(PAMcand)
    if not d1.empty:
    # Check if the target sequence starts with the spacer
        if target[:len(spacer)] == spacer:
            pamDir = 3
            d1['pam'] = d1['Site_SubstitutionsOnly.Sequence'].apply(lambda x: x[len(spacer):len(target)])  # 3' PAM
        else:
            pamDir = 5
            d1['pam'] = d1['Site_SubstitutionsOnly.Sequence'].apply(lambda x: x[:len(target)-len(spacer)])  # 5' PAM
            
            # Print PAM Direction
        print("PAM Direction:", pamDir)
        print(d1)
        print(d1['pam'])
    d2 = d1[d1['pam'] != ""].copy()

    # Assign values to column 'PMMM' based on conditions
    d2['PMMM'] = np.where(d2['Site_SubstitutionsOnly.NumSubstitutions'] == 0, 'PM', 
                np.where((d2['Site_SubstitutionsOnly.NumSubstitutions'] > 0) & 
                (d2['Site_SubstitutionsOnly.NumSubstitutions'] < 7), 'MM', ''))
    print(d2)
        # Loop through PM and MM types
    for mType in ['PM', 'MM']:
        dataPMMM = d2[d2['PMMM'] == mType]
        nsites = len(dataPMMM)
        if nsites != 0:
            pos_start = 1
            pos_pam = np.nan
            pam_d = np.nan
            peak = 1
            for i in range(lenPAM - 3):
                pam_i = np.repeat(dataPMMM['pam'].str[i:i+3], dataPMMM['bi.sum.mi'])
                pam_tabi = pd.Series(pam_i).value_counts()
                pam_maxi = pam_tabi.max()
                sum_im = pam_tabi.sum()
                print(f"i: {i}; max cnt: {pam_maxi}; sum: {sum_im}; pam_table: {pam_tabi}")
                if pam_maxi > peak:
                    peak = pam_maxi
                    pos_start = i
                    pam_d = pd.DataFrame({'pam_i': pam_tabi.index})
                    pam_d['x1'] = pam_d['pam_i'].str[0]
                    pam_d['x2'] = pam_d['pam_i'].str[1]
                    pam_d['x3'] = pam_d['pam_i'].str[2]
                    pam_d['x4'] = pam_d['pam_i'].str[3]

                    for x1 in ['A', 'G', 'C', 'T']:
                        for x2 in ['A', 'G', 'C', 'T']:
                            for x3 in ['A', 'G', 'C', 'T']:
                                for x4 in ['A', 'G', 'C', 'T']:
                                    pam_d = pam_d.append(pd.Series([f"{x1}{x2}{x3}{x4}", 0, x1, x2, x3, x4], index=pam_d.columns), ignore_index=True)
                    pam_d = pam_d[~pam_d['pam_i'].duplicated()]
                    pam_d['Fraction'] = pd.to_numeric(pam_d['Freq']) / pam_d['Freq'].sum()
                    pos_str = pos_start if pamDir == 5 else pos_start - lenPAM - 1
                    pam_d['x1'] = pam_d['x1'].apply(lambda x: f"Position {pos_str}: {x}")
                    pam_d['x2'] = pam_d['x2'].apply(lambda x: f"Position {pos_str + 1}: {x}")

                    # Insert the snippet here
                    # Filter out empty strings in pam_d DataFrame
                    pam_d_filtered = pam_d[pam_d['pam_i'] != '']
                    pam_d_filtered['pam_i'] = pam_d_filtered['pam_i'].apply(lambda x: x[0] if len(x) > 0 else '')
                    print(pam_d_filtered)
                    plt.figure(figsize=(10, 8))

                    # Create the plot
                    sns.heatmap(data=pam_d.pivot(index='x2', columns='x1', values='Fraction'), cmap='coolwarm', annot=False)

                    # Set the title and axis labels
                    plt.title(f"{runDate}_{LibID}_PAM\n({nsites} {mType} sites) positions: {pos_str} to {pos_str+3}")
                    plt.xlabel(f"Position {pos_str+2}")
                    plt.ylabel(f"Position {pos_str+3}")

                    # Save the plot as a PDF
                    plt.tight_layout()
                    plt.savefig(f"{runDate}_{LibID}_PAM_{mType}_{pos_str}-{pos_str+3}.pdf")

                    # Display the plot
                    plt.show()
def main():
    IdentifiedOfftargetsFile = sys.argv[1]
    runDate = sys.argv[2]
    generate_dataframe(IdentifiedOfftargetsFile,runDate)

if __name__ == "__main__":
    main()

# def plot_nice_heatmap(df_output, avg_spacer, spacer=None):
#     ## Heatmap plotting
#     axes = [4 * (pam_length - split_pam_index), 4 * (split_pam_index)]
#     fig, ax1 = plt.subplots(1, figsize=(axes[0], axes[1]))
#     sns.heatmap(df_output,
#                 vmin=0,
#                 vmax=1,
#                 square=False,
#                 cman=''