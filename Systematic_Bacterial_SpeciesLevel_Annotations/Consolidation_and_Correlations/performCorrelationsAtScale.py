import os
import sys
from scipy import stats
from collections import defaultdict
import pandas as pd

data_file = 'DataFrame.txt'
df = pd.read_csv(data_file, sep='\t', header=0)

#y_train = df2['bgc_sum']

bac_df = df_full.loc[df_full['phylum'] == 'Bacillota']
act_df = df_full.loc[df_full['phylum'] == 'Actinomycetota']
pse_df = df_full.loc[df_full['phylum'] == 'Pseudomonadota']
myx_df = df_full.loc[df_full['phylum'] == 'Myxococcota']
cya_df = df_full.loc[df_full['phylum'] == 'Cyanobacteriota']

dfs = [df_full, bac_df, act_df, pse_df, myx_df, cya_df]
names = ['General', 'Bacillota', 'Actinomycetota', 'Pseudomonadota', 'Myxococcota', 'Cyanobacteriota']

aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

for i, df in enumerate(dfs):
    df_name = names[i]

    bgc_sums = df['bgc_sum']
    genome_sizes = df['genome_size']
    r2, pval = stats.spearmanr(bgc_sums, genome_sizes)
    print(df_name + '\tGenome Size\t' + str(r2) + '\t' + str(pval))  

    transposon_counts = df['transposon_count']
    r2, pval = stats.spearmanr(bgc_sums, transposon_counts)
    print(df_name + '\tTransposon Homologs\t' + str(r2) + '\t' + str(pval))

    cazy_counts = df['cazy_count'] 
    r2, pval = stats.spearmanr(bgc_sums, cazy_counts)
    print(df_name + '\tDistinct CAZy Enzyme Homologs\t' + str(r2) + '\t' + str(pval))

    oxy_counts = df['oxidative_phospho_count']
    r2, pval = stats.spearmanr(bgc_sums, oxy_counts)
    print(df_name + '\tDistinct Oxidative Phosphorylation Homologs\t' + str(r2) + '\t' + str(pval))

    adh_counts = df['adherence_count']
    r2, pval = stats.spearmanr(bgc_sums, adh_counts)
    print(df_name + '\tDistinct Adherence Protein Homologs\t' + str(r2) + '\t' + str(pval))

    gc = df['gc']
    r2, pval = stats.spearmanr(bgc_sums, gc)
    print(df_name + '\tGC%\t' + str(r2) + '\t' + str(pval))

    ideal_temps = df['ideal_temp']
    r2, pval = stats.spearmanr(bgc_sums, ideal_temps)
    print(df_name + '\tInferred Ideal Growth Temperature\t' + str(r2) + '\t' + str(pval))

    phage_sums = df['phage_sum']
    r2, pval = stats.spearmanr(bgc_sums, phage_sums)
    print(df_name + '\tPhage-ome Size\t' + str(r2) + '\t' + str(pval))

    bgc_sums = df.loc[df['host_associated_category'] != 'not-detected']['bgc_sum']
    med_rel_abds = df.loc[df['host_associated_category'] != 'not-detected']['median_relative_abundance']
    r2, pval = stats.spearmanr(bgc_sums, med_rel_abds)
    print(df_name + '\tMedian Relative Abundance in Metagenomes\t' + str(r2) + '\t' + str(pval))

    prop_host_associated = df.loc[df['host_associated_category'] != 'not-detected']['proportion_host_associated']
    r2, pval = stats.spearmanr(bgc_sums, prop_host_associated)
    print(df_name + '\tProportion of Metagenomes Detected Eukaryotic-Host Associated\t' + str(r2) + '\t' + str(pval))

    host_mgs = df.loc[df['host_associated_category'] != 'not-detected']['host_metagenome_count']
    r2, pval = stats.spearmanr(bgc_sums, host_mgs)
    print(df_name + '\tHost Metagenomes Detected Within\t' + str(r2) + '\t' + str(pval))

    env_mgs = df.loc[df['host_associated_category'] != 'not-detected']['environment_metagenome_count']
    r2, pval = stats.spearmanr(bgc_sums, env_mgs)
    print(df_name + '\tEnvironment Metagenomes Detected Within\t' + str(r2) + '\t' + str(pval))

    tot_mgs = df.loc[df['host_associated_category'] != 'not-detected']['total_metagenome_count']
    r2, pval = stats.spearmanr(bgc_sums, tot_mgs)
    print(df_name + '\tTotal Metagenomes Detected Within\t' + str(r2) + '\t' + str(pval))
    
    bgc_sums = df['bgc_sum']
    for aa in aminos:
        aa_freqs = df[aa]
        r2, pval = stats.spearmanr(bgc_sums, aa_freqs)
        print(df_name + '\t' + aa + ' Amino Acid Frequency\t' + str(r2) + '\t' + str(pval))
    
