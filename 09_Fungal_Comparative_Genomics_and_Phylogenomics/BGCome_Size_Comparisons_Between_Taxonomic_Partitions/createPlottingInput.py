import os
import sys
from scipy import stats
from collections import defaultdict

antismash_stats_file = 'AntiSMASH_Stats.Updated_WithoutContaminants.txt'

pez_clade_file = 'Clades/Pezizomycotina.txt'
yea_clade_file = 'Clades/Dikarya_Yeasts.txt'
bas_clade_file = 'Clades/NonYeast_Basidi.txt'
neo_clade_file = 'Clades/Neocalli.txt'
zoo_clade_file = 'Clades/Zoopago.txt'
muc_clade_file = 'Clades/Mucoro.txt'
spi_clade_file = 'Clades/Spizello.txt'
out_clade_file = 'Clades/Other_Zoosporic_or_Basal.txt'

genus_reps_file = 'Genus_Representative_GCAs.txt'
genus_reps = set([])
with open(genus_reps_file) as ogrf:
    for line in ogrf:
        line = line.strip()
        genus_reps.add(line)


clade_lists = [pez_clade_file, bas_clade_file, yea_clade_file, muc_clade_file, zoo_clade_file, neo_clade_file, spi_clade_file, out_clade_file]
clade_names = ['1. Pezizomycotina', '2. Non-Yeast Basidimycota', '3. Yeast Dikarya', '4. Mucoromycota', '5. Zoopagomycota', '6. Neocallimastigomycota', '7. Spizellomycetales', '8. Other Zoosporic and Basal Fungi']
colors = ['haploid-dominant', 'polyploid-dominant', 'yeast', 'haploid-dominant', 'polyploid-dominant', 'haploid-dominant', 'haploid-dominant', 'polyploid-dominant']

clade_gcas = defaultdict(set)
gca_to_clade = {}
gca_to_color = {}
for i, clf in enumerate(clade_lists):
    cn = clade_names[i]
    color = colors[i]
    with open(clf) as ocf:
        for line in ocf:
            line = line.strip()
            gca = '_'.join(line.split('_')[-2:])
            if gca in genus_reps:
                clade_gcas[cn].add(gca)
                gca_to_clade[gca] = cn
                gca_to_color[gca] = color

print('GCA\tPhylum_or_Clade\tComplete_BGC_Count\tBGCome_Size\tStrict_NRPS_or_PKSome_Size\tColor\tGrid')
clade_bgcome_sizes = defaultdict(list)
with open(antismash_stats_file) as oasf:
    for line in oasf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        if not gca in gca_to_clade: continue
        clade = gca_to_clade[gca]
        color = gca_to_color[gca]
        bgcome_size = float(ls[4])/1000000.0
        comp_bgc_count = float(ls[3])
        strict_nrps_or_pksome_size = float(ls[-1])/1000000.0
        clade_bgcome_sizes[clade].append(comp_bgc_count)#strict_nrps_or_pksome_size*1000000.0)
        grid = '2 other'
        if clade == '1. Pezizomycotina':
            grid = '1 Pez'
        print('\t'.join([gca, clade, str(comp_bgc_count), str(bgcome_size),  str(strict_nrps_or_pksome_size), color, grid]))

print('Non-Yeast Basidi vs. All Yeast Dikarya')
print(stats.ranksums(clade_bgcome_sizes['2. Non-Yeast Basidimycota'], clade_bgcome_sizes['3. Yeast Dikarya'], alternative='greater'))

print('Muco vs. Zoop')
print(stats.ranksums(clade_bgcome_sizes['4. Mucoromycota'], clade_bgcome_sizes['5. Zoopagomycota'], alternative='greater'))

print('Neo vs. Other Zoosporic + Basal')
print(stats.ranksums(clade_bgcome_sizes['6. Neocallimastigomycota'], clade_bgcome_sizes['8. Other Zoosporic and Basal Fungi'], alternative='greater'))

print('Spiz & Rhiz vs. Other Zoosporic + Basal')
print(stats.ranksums(clade_bgcome_sizes['7. Spizellomycetales'], clade_bgcome_sizes['8. Other Zoosporic and Basal Fungi'], alternative='greater'))

