import os
import sys
from collections import defaultdict

gca_to_name = {}
with open('../AntiSMASH_Stats.Updated_WithoutContaminants.txt') as ofgg:
    for line in ofgg:
        line = line.strip()
        ls = line.split('\t')
        gca_to_name['_'.join(ls[1].split('_')[-2:])] = ls[1]

og_hits = defaultdict(lambda: defaultdict(int))
gca_hits = defaultdict(int)

assoc_ogs = set([])
with open('HetAssocOGs.txt') as ohf:
    for line in ohf:
        line = line.strip()
        ls = line.split('\t')
        assoc_ogs.add(ls[0])

gcas = []
with open(sys.argv[1]) as of:
    for i, line in enumerate(of):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            gcas = ls[1:]
            continue
        og = ls[0]
        if not og in assoc_ogs: continue
        for j, cdss in enumerate(ls[1:]):
            gca = gcas[j]
            for cds in cdss.split(', '):
                if cds.strip() != '':
                    gca_hits[gca] += 1

print('DATASET_SIMPLEBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tHIP')
print('COLOR\t#000000')
print('DATA')
for gca in gca_hits:
    print(gca_to_name[gca] + '\t' + str(gca_hits[gca]))

