import os
import sys
from collections import defaultdict
from operator import itemgetter

og_file = 'Orthogroups.tsv'
map_file = '../Representative_Genomes_for_Comparative_Genomics.txt'
sog_file = 'Select_OGs.txt'

name_map = {}
with open(map_file) as omf:
    for line in omf:
        line = line.strip()
        ls = line.split('\t')
        name_map[ls[-1].split('.')[0]] = ls[0]

og_bgc_genera = defaultdict(set)
with open(sog_file) as osf:
    for line in osf:
        line = line.strip('\n')
        ls = line.split('\t')
        for g in ls[3].split(', '):
            og_bgc_genera[ls[0]].add(g)
        for g in ls[4].split(', '):
            og_bgc_genera[ls[0]].add(g)
        for g in ls[5].split(', '):
            og_bgc_genera[ls[0]].add(g)

genera = []
genus_sample_ogs = defaultdict(lambda: defaultdict(int))
og_counts = defaultdict(int)
with open(og_file) as oof:
    for i, line in enumerate(oof):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            genera = [name_map[x.split('.')[0]] for x in ls[1:]]
        else:
            og = ls[0]
            if og in og_bgc_genera:
                for j, val in enumerate(ls[1:]):
                    genus = genera[j]
                    if val == '': continue
                    if genus in og_bgc_genera[og]:
                        genus_sample_ogs[og][genus] = 2
                    else:
                        genus_sample_ogs[og][genus] = 1
                    og_counts[og] += 1

og_order = []
for og in sorted(og_counts.items(), key=itemgetter(1), reverse=True):
   og_order.append(og[0])

print('DATASET_HEATMAP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tBGC Key OGs')
print('FIELD_LABELS\t' + '\t'.join(og_order))
print('DATA')
for g in genera:
    printlist = [g]
    for og in og_order:
        printlist.append(str(genus_sample_ogs[og][g]))
    print('\t'.join(printlist))

