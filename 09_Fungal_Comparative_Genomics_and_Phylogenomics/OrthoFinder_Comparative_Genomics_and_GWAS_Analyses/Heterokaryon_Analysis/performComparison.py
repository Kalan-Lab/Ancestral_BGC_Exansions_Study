import os
import sys
from collections import defaultdict
from scipy import stats

bgc_enriched_clade = '/home/salamzade/Multicellularity_and_BGC_Abundances/Polished_Fungi/antiSMASH/BGCome_Sizes/Clades/BGC_Enriched_Clade.txt'
dikarya_clade = '/home/salamzade/Multicellularity_and_BGC_Abundances/Polished_Fungi/antiSMASH/BGCome_Sizes/Clades/Dikarya.txt'


bgc_enriched = set([])
dikarya = set([])

with open(bgc_enriched_clade) as obec:
    for line in obec:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        bgc_enriched.add(gca)

with open(dikarya_clade) as odc:
    for line in odc:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        dikarya.add(gca)

og_hits = defaultdict(lambda: defaultdict(int))
gca_hits = defaultdict(int)

hit_cds = set([])
with open('Heterokaryon_Search_All_Prots.txt') as of:
    for line in of:
        line = line.strip()
        if line.startswith('#'): continue
        ls = line.split()
        evalue = float(ls[4])
        if evalue < 1e-3 and not ls[0] in hit_cds:
            hit_cds.add(ls[0])
            gca_hits[ls[0].split('|')[0]] += 1

bec = []
odc = []
for gca in gca_hits:
    if gca in dikarya:
        if gca in bgc_enriched:
            bec.append(gca_hits[gca])
        else:
            odc.append(gca_hits[gca])

print(stats.ranksums(bec, odc, alternative='greater'))
