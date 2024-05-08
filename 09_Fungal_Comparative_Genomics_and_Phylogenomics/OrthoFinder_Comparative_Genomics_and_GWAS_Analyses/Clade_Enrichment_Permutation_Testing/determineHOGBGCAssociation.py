import os
import sys
from collections import defaultdict

bgc_file = 'All_BGC_Proteins.with_Labels.txt'
hog_file = 'Hierarchical_Orthogroups.tsv'

hog_bgc_counts = defaultdict(int)
hog_key_bgc_counts = defaultdict(int)
bgc_region_found = set([])
key_bgc_hog = set([])
with open(bgc_file) as obf:
    for line in obf:
        line = line.strip()
        ls = line.split('\t')
        if not 'N0' in ls[0]: continue
        hog = ls[0].split('N0.')[1]
        hog_bgc_counts[hog] += 1
        bgc_region_found.add(hog)
        if ls[-1] == 'True':
            key_bgc_hog.add(hog)
            hog_key_bgc_counts[hog] += 1
    
total_hog_counts = defaultdict(int)
with open(hog_file) as ohf:
    for i, line in enumerate(ohf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: continue
        hog = ls[0].split('N0.')[1]
        if hog in bgc_region_found:
            for gs in ls[1:]:
                for g in gs.split(', '):
                    if g.strip() != '':
                        total_hog_counts[hog] += 1

for hog in bgc_region_found:
    print(hog + '\t' + str(total_hog_counts[hog]) + '\t' + str(hog_bgc_counts[hog]/total_hog_counts[hog]) + '\t' + str(hog_key_bgc_counts[hog]/total_hog_counts[hog]))
