import os
import sys
from collections import defaultdict
from scipy import stats

bep_clade_file = 'Clades/BGC_Enriched_Pezizomycotina.txt'
oth_clade_file = 'Clades/Non_BGCrichPezi_Fungi.txt'
ndi_clade_file = 'Clades/NonDikarya.txt'
pez_clade_file = 'Clades/Pezizomycotina.txt'
asc_clade_file = 'Clades/Ascomycota_Yeasts_and_Other.txt'
bas_clade_file = 'Clades/Basidiomycota_Yeasts_and_Other.txt'
aga_clade_file = 'Clades/Agaricomycetes.txt'
neo_clade_file = 'Clades/Neocallimastigomycota.txt'
zoo_clade_file = 'Clades/Zoopagomycota.txt'
muc_clade_file = 'Clades/Mucoromycota.txt'
chy_clade_file = 'Clades/Chytrid_and_ChytridLike.txt'
bla_clade_file = 'Clades/Blastocladiomycota.txt'
out_clade_file = 'Clades/Basal_Fungi.txt'
chl_clade_file = 'Clades/ChytridAndRelated.txt'
sac_clade_file = 'Clades/Saccharomycotina.txt'
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

clade_lists = [bep_clade_file, pez_clade_file, asc_clade_file, aga_clade_file, bas_clade_file, muc_clade_file, zoo_clade_file, neo_clade_file, chy_clade_file, bla_clade_file, out_clade_file, oth_clade_file, ndi_clade_file, chl_clade_file, sac_clade_file, nco_clade_file]
clade_names = ['BGC-enriched Pezizomycotina', 'Pezizomycotina', 'Yeast and Other Ascomycota', 'Agaricomycetes', 'Other Basidiomycota', 'Mucoromycota', 'Zoopagomycota', 'Neocallimastigomycota', 'Chytrid and Other Chytrid-Like', 'Blastocladiomycota', 'Basal Fungi', 'Other fungi', 'non-Dikaryon fungi', 'Chytrid & Chytrid-like (including Neocallimastigomycota)', 'Saccharomycotina', 'Other non-Dikaryon fungi']

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

clade_comparisons = [['BGC-enriched Pezizomycotina', 'Other fungi'], ['Pezizomycotina', 'Agaricomycetes'], ['Agaricomycetes', 'Other Basidiomycota'], ['Agaricomycetes', 'non-Dikaryon fungi'], ['Chytrid & Chytrid-like (including Neocallimastigomycota)', 'Saccharomycotina'], ['Neocallimastigomycota', 'Other non-Dikaryon fungi']]

outf = open('Plotting_Input.txt', 'w')
outf.write('\t'.join(['clade_comparison', 'clade', 'genome', 'bgcome_size']) + '\n')
for i, cc in enumerate(clade_comparisons):

    cc_name = str(i) + '. ' + ' vs. '.join(cc)
    c1, c2 = cc
    c1_bs = clade_bgcome_sizes[c1]
    c2_bs = clade_bgcome_sizes[c2]
    print(cc_name)
    print(str(len(c1_bs)) + '\t' + str(len(c2_bs)))
    stat, pval = stats.ranksums(c1_bs, c2_bs, alternative='greater')
    print(pval)
    print('-------------------')
    for g in clade_gcas[c1]:
        gbs = bgcome_sizes[g]
        outf.write('\t'.join([cc_name, '1. ' + c1, g, str(gbs)]) + '\n')
    for g in clade_gcas[c2]:
        gbs = bgcome_sizes[g]
        outf.write('\t'.join([cc_name, '2. ' + c2, g, str(gbs)]) + '\n')

outf.close()
