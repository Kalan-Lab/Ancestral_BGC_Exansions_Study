import os
import sys
from collections import defaultdict

pks_annot_file = 'PKS_Annotation.txt'
rifk_to_gc_file = '../RifK_to_FaiGC.txt'

mapping = {}
with open(rifk_to_gc_file) as orgf:
    for line in orgf:
        line = line.strip()
        ls = line.split('\t')
        mapping[ls[1].replace('.gbk', '')] = ls[0]

gc_pks_sets = defaultdict(set)
with open(pks_annot_file) as opaf:
    for line in opaf:
        line = line.strip()
        ls = line.split('\t')
        slen = int(ls[-1])
        gc = ls[1].split('|')[0]
        if slen >= 1000:
            gc_pks_sets[gc].add(ls[1])


print('DATASET_SIMPLEBAR')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tPKS Counts')
print('DATA')
for gc in gc_pks_sets:
    if gc in mapping:
        print(mapping[gc] + '\t' + str(len(gc_pks_sets[gc])))
    #print(gc + '\t' + str(len(gc_pks_sets[gc])))
