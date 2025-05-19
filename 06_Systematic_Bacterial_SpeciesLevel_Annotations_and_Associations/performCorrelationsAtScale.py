import os
import sys
from scipy import stats
from collections import defaultdict
import pandas as pd

data_file = 'Comprehensive_Info_for_Bacterial_Species_Reps.tsv'
df_full = pd.read_csv(data_file, sep='\t', header=0)

bac_df = df_full.loc[df_full['phylum'] == 'Bacillota']
act_df = df_full.loc[df_full['phylum'] == 'Actinomycetota']
myx_df = df_full.loc[df_full['phylum'] == 'Myxococcota']
cya_df = df_full.loc[df_full['phylum'] == 'Cyanobacteriota']
pse_df = df_full.loc[df_full['phylum'] == 'Pseudomonadota']

dfs = [df_full, pse_df, bac_df, cya_df, myx_df, act_df]
names = ['General', 'Pseudomonadota', 'Bacillota', 'Cyanobacteriota', 'Myxococcota', 'Actinomycetota']
 
aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

main_metric = 'strict_nrps_or_pks_ome_size' # 'relaxed_bgcome_size'

for i, df in enumerate(dfs):
    df_name = names[i]

    y = df[main_metric]

    genome_sizes = df['genome_size']
    r2, pval = stats.spearmanr(y, genome_sizes)
    print(df_name + '\tGenome Size\t' + str(r2) + '\t' + str(pval))  

    transposon_counts = df['transposon_count']
    r2, pval = stats.spearmanr(y, transposon_counts)
    print(df_name + '\tTransposon Homologs\t' + str(r2) + '\t' + str(pval))

    cazy_distinct_counts = df['cazy_distinct_count'] 
    r2, pval = stats.spearmanr(y, cazy_distinct_counts)
    print(df_name + '\tDistinct CAZy Enzyme Families\t' + str(r2) + '\t' + str(pval))

    cazy_total_counts = df['cazy_total_count']
    r2, pval = stats.spearmanr(y, cazy_total_counts)
    print(df_name + '\tTotal CAZy Enzyme Homologs\t' + str(r2) + '\t' + str(pval))

    oxy_counts = df['oxidative_phospho_count']
    r2, pval = stats.spearmanr(y, oxy_counts)
    print(df_name + '\tDistinct Oxidative Phosphorylation Homologs\t' + str(r2) + '\t' + str(pval))

    adh_counts = df['adherence_count']
    r2, pval = stats.spearmanr(y, adh_counts)
    print(df_name + '\tDistinct Adherence Protein Homologs\t' + str(r2) + '\t' + str(pval))

    gc = df['gc']
    r2, pval = stats.spearmanr(y, gc)
    print(df_name + '\tGC%\t' + str(r2) + '\t' + str(pval))

    ideal_temps = df['ideal_temp']
    r2, pval = stats.spearmanr(y, ideal_temps)
    print(df_name + '\tPredicted Optimal Growth Temperature\t' + str(r2) + '\t' + str(pval))

    phage_sums = df['phage_sum']
    r2, pval = stats.spearmanr(y, phage_sums)
    print(df_name + '\tPhage-ome Size\t' + str(r2) + '\t' + str(pval))

    y = df.loc[df['host_associated_category'] != 'not-detected'][main_metric]
    
    med_rel_abds = df.loc[df['host_associated_category'] != 'not-detected']['median_relative_abundance']
    r2, pval = stats.spearmanr(y, med_rel_abds)
    print(df_name + '\tMedian Relative Abundance in MGs\t' + str(r2) + '\t' + str(pval))

    prop_host_associated = df.loc[df['host_associated_category'] != 'not-detected']['proportion_host_associated']
    r2, pval = stats.spearmanr(y, prop_host_associated)
    print(df_name + '\tProp. of MGs which are Eukaryotic-Host Associated\t' + str(r2) + '\t' + str(pval))

    #host_mgs = df.loc[df['host_associated_category'] != 'not-detected']['host_metagenome_count']
    #r2, pval = stats.spearmanr(bgc_sums, host_mgs)
    #print(df_name + '\tHost MGs Detected Within\t' + str(r2) + '\t' + str(pval))

    #env_mgs = df.loc[df['host_associated_category'] != 'not-detected']['environment_metagenome_count']
    #r2, pval = stats.spearmanr(bgc_sums, env_mgs)
    #print(df_name + '\tEnvironment Metagenomes Detected Within\t' + str(r2) + '\t' + str(pval))

    tot_mgs = df.loc[df['host_associated_category'] != 'not-detected']['total_metagenome_count']
    r2, pval = stats.spearmanr(y, tot_mgs)
    print(df_name + '\tTotal MGs Detected Within\t' + str(r2) + '\t' + str(pval))
    
    y = df[main_metric]

    for aa in aminos:
        aa_freqs = df[aa]
        r2, pval = stats.spearmanr(y, aa_freqs)
        print(df_name + '\t' + aa + ' Amino Acid Frequency\t' + str(r2) + '\t' + str(pval))
    
