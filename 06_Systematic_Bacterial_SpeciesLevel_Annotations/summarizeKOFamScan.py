import os
import sys
from collections import defaultdict

results_dir = 'KOFamScan_Select_Results/'
spl_sel_file = 'Species_Level_Selections.txt'

ko_nf = set(['K02588', 'K02586', 'K02591', 'K00531', 'K22896', 'K22897', 'K22898', 'K22899'])

gca_to_taxa = {}
with open(spl_sel_file) as ossf:
    for line in ossf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[-1].split('.')[0]
        gca_to_taxa[gca] = ls[0] + '\t' + ls[1]

print('GCA\tPhylum\tGenus\tKO\tKO_Count\tKO_Category')
for f in os.listdir(results_dir):
    gca = f.split('.')[0]
    ko_count = defaultdict(int)
    with open(results_dir + f) as orf:
        for line in orf:
            line = line.strip()
            ls = line.split()
            if line.startswith('*'): 
                ko_count[ls[2]] += 1
    for ko in ko_count:
        ko_cat = 'Oxidative Phosphorylation'
        if ko in ko_nf:
            ko_nf = 'Nitrogen Fixation'
        print(gca + '\t' + gca_to_taxa[gca] + '\t' + ko + '\t' + str(ko_count[ko]) + '\t' + ko_cat)
