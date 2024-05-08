import os
import sys
from collections import defaultdict

gca_to_id = {}
with open('../AntiSMASH_Stats.Updated_WithoutContaminants.txt') as oaf:
    for i, line in enumerate(oaf):
        line = line.strip()
        ls = line.split('\t')
        gca_to_id[ls[0]] = ls[1]

select_ogs = set([])
with open('Yeast_Prions_Mapping.txt') as ogmf:
    for line in ogmf:
        line = line.strip()
        ls = line.split('\t')
        select_ogs.add(ls[-1])

gca_og_counts = defaultdict(int)
gcas = []
with open('Orthogroups.tsv') as ogf:
    for i, line in enumerate(ogf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            gcas = ls[1:]
        else:
            og = ls[0]
            if not og in select_ogs: continue
            for j, gs in enumerate(ls[1:]):
                gca = gcas[j]
                if gs.strip() != '':
                    for g in gs.split(', '):
                        gca_og_counts[gca] += 1 # len(gs.split(', '))


print('DATASET_SIMPLEBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tYeast prion homologs')
print('COLOR\t#000000')
print('DATA')
for gca in gcas:
    print(gca_to_id[gca] + '\t' + str(gca_og_counts[gca]))
