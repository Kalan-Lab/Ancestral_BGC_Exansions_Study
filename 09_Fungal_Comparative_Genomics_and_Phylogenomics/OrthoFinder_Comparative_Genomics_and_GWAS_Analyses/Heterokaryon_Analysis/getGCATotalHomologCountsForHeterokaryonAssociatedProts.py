import os
import sys
from collections import defaultdict

gca_to_name = {}
with open('AntiSMASH_Stats.Updated_WithoutContaminants.txt') as ofgg:
    for line in ofgg:
        line = line.strip()
        ls = line.split('\t')
        gca_to_name['_'.join(ls[1].split('_')[-2:])] = ls[1]

og_hits = defaultdict(lambda: defaultdict(int))
gca_hits = defaultdict(int)

hit_cds = set([])
with open('Heterokaryon_Search_All_Prots.txt') as of:
    for line in of:
        line = line.strip()
        if line.startswith('#'): continue
        ls = line.split()
        evalue = float(ls[4])
        if evalue < 1e-3:
            hit_cds.add(ls[0])
            gca_hits[ls[0].split('|')[0]] += 1

"""
gcas = []
with open(sys.argv[1]) as of:
    for i, line in enumerate(of):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            gcas = ls[1:]
            continue
        og = ls[0]
        for j, cdss in enumerate(ls[1:]):
            gca = gcas[j]
            for cds in cdss.split(', '):
                if cds.strip() != '' and cds in hit_cds:
                    og_hits[og][gca] += 1
                    gca_hits[gca] += 1
"""

print('DATASET_SIMPLEBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tHIP')
print('COLOR\t#000000')
print('DATA')

for gca in gca_hits:
    print(gca_to_name[gca] + '\t' + str(gca_hits[gca]))

#for og in og_hits:
#    for gca in og_hits[og]:
#        print(og + '\t' + gca + '\t' + str(og_hits[og][gca]))
