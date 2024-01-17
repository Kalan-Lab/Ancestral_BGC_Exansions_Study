import os
import sys
from scipy import stats
from collections import defaultdict
import pandas as pd

"""
0       gca
1       clade
2       bgc_sum
3       bgc_prop
4       genome_size
5       gc
6       cazy_counts
7       starship_counts
8       A
9       C
10      D
11      E
12      F
13      G
14      H
15      I
16      K
17      L
18      M
19      N
20      P
21      Q
22      R
23      S
24      T
25      V
26      W
27      Y
"""


data_file = 'DataFrame.txt'
df_full = pd.read_csv(data_file, sep='\t', header=0)

enr_df = df_full.loc[df_full['clade'] == 'Ascomycota - BGC Enriched Clade']

dfs = [df_full, enr_df]
names = ['General', 'Ascomycota - BGC Enriched Clade']

aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

for i, df in enumerate(dfs):
    df_name = names[i]

    bgc_sums = df['bgc_sum']
    genome_sizes = df['genome_size']
    r2, pval = stats.spearmanr(bgc_sums, genome_sizes)
    print(df_name + '\tGenome Size\t' + str(r2) + '\t' + str(pval))  

    transposon_counts = df['starship_count']
    r2, pval = stats.spearmanr(bgc_sums, transposon_counts)
    print(df_name + '\tStarship Conserved Gene Homologs\t' + str(r2) + '\t' + str(pval))

    cazy_counts = df['cazy_count']
    r2, pval = stats.spearmanr(bgc_sums, cazy_counts)
    print(df_name + '\tDistinct CAZy Enzyme Homologs\t' + str(r2) + '\t' + str(pval))

    gc = df['gc']
    r2, pval = stats.spearmanr(bgc_sums, gc)
    print(df_name + '\tGC%\t' + str(r2) + '\t' + str(pval))

    bgc_sums = df['bgc_sum']
    for aa in aminos:
        aa_freqs = df[aa]
        r2, pval = stats.spearmanr(bgc_sums, aa_freqs)
        print(df_name + '\t' + aa + ' Amino Acid Frequency\t' + str(r2) + '\t' + str(pval))
    
