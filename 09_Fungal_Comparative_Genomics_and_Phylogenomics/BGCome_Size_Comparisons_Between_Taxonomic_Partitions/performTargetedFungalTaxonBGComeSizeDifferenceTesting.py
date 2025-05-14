import os
import sys
from collections import defaultdict
from scipy import stats

bep_clade_file = 'Clades/BGC_Enriched_Pezizomycotina.txt'
oth_clade_file = 'Clades/Non_BGCrichPezi_Fungi.txt'
ndi_clade_file = 'Clades/NonDikarya.txt'
bas_clade_file = 'Clades/Basidiomycota_Yeasts_and_Other.txt'
aga_clade_file = 'Clades/Agaricomycetes.txt'
neo_clade_file = 'Clades/Neocallimastigomycota.txt'
nco_clade_file = 'Clades/NonDikaryonNonNeo.txt'
genus_reps_file = 'Representatives.txt'

genus_reps = set([])
with open(genus_reps_file) as ogrf:
    for line in ogrf:
        line = line.strip()
        genus_reps.add(line.split('\t')[0])

antismash_stats_file = 'AntiSMASH_Stats.txt'
bgcome_sizes = {}
with open(antismash_stats_file) as oasf:
    for i, line in enumerate(oasf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        name = ls[1]
        bgcome_sizes[name] = float(ls[4])

clade_lists = [bep_clade_file, oth_clade_file, ndi_clade_file, bas_clade_file, aga_clade_file, neo_clade_file, nco_clade_file]
clade_names = ['1. BGC-enriched Pezizomycotina', '2. Other fungi', '3. Non-Dikaryon fungi', '2. Other Basidiomycota', '1. Agaricomycetes', '1. Neocallimastigomycota', '2. Other non-Dikaryon fungi']

clade_gcas = defaultdict(set)
clade_bgcome_sizes = defaultdict(list)
for i, clf in enumerate(clade_lists):
    cn = clade_names[i]
    with open(clf) as ocf:
        for line in ocf:
            line = line.strip()
            name = line
            gca = '_'.join(line.split('_')[1:]).replace('GCF_', 'GCA_').split('.')[0]
            if gca in genus_reps: 
                clade_gcas[cn].add(name)
                clade_bgcome_sizes[cn].append(bgcome_sizes[name])

clade_comparisons = [['1. BGC-enriched Pezizomycotina', '2. Other fungi'], ['1. Neocallimastigomycota', '2. Other non-Dikaryon fungi'], ['1. Agaricomycetes', '2. Other Basidiomycota', '3. Non-Dikaryon fungi']]

outf = open('Plotting_Input.txt', 'w')
outf.write('\t'.join(['clade_comparison', 'clade', 'genome', 'bgcome_size']) + '\n')
for i, cc in enumerate(clade_comparisons[:-1]):

    cc_name = str(i) + '. ' + ' vs. '.join(cc)
    c1, c2 = cc
    c1_bs = clade_bgcome_sizes[c1]
    c2_bs = clade_bgcome_sizes[c2]
    #print(cc_name)
    #print(str(len(c1_bs)) + '\t' + str(len(c2_bs)))
    stat, pval = stats.ranksums(c1_bs, c2_bs, alternative='greater')
    #print(pval)
    #print('-------------------')
    for g in clade_gcas[c1]:
        gbs = bgcome_sizes[g]
        outf.write('\t'.join([cc_name, c1, g, str(gbs)]) + '\n')
    for g in clade_gcas[c2]:
        gbs = bgcome_sizes[g]
        outf.write('\t'.join([cc_name, c2, g, str(gbs)]) + '\n')

i = 2
c1, c2, c3 = clade_comparisons[-1]
c1_bs = clade_bgcome_sizes[c1]
c2_bs = clade_bgcome_sizes[c2]
c3_bs = clade_bgcome_sizes[c3]

#print('Agari comps')
stat, pval = stats.ranksums(c1_bs, c2_bs, alternative='greater')
#print(pval)
stat, pval = stats.ranksums(c1_bs, c3_bs, alternative='greater')
#print(pval)
#print('-------------------')

for g in clade_gcas[c1]:
    gbs = bgcome_sizes[g]
    outf.write('\t'.join(['Agari comps', c1, g, str(gbs)]) + '\n')
for g in clade_gcas[c2]:
    gbs = bgcome_sizes[g]
    outf.write('\t'.join(['Agari comps', c2, g, str(gbs)]) + '\n')
for g in clade_gcas[c3]:
    gbs = bgcome_sizes[g]
    outf.write('\t'.join(['Agari comps', c3, g, str(gbs)]) + '\n')
outf.close()
