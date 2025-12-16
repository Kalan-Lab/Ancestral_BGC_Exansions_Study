import os
import sys
from collections import defaultdict

is_annot_file = 'IS_Annotation.txt'
rifk_to_gc_file = '../RifK_to_FaiGC.txt'

mapping = {}
with open(rifk_to_gc_file) as orgf:
    for line in orgf:
        line = line.strip()
        ls = line.split('\t')
        mapping[ls[1].split('.gbk')[0]] = ls[0]

gc_is_sets = defaultdict(set)
with open(is_annot_file) as opaf:
    for line in opaf:
        line = line.strip()
        ls = line.split('\t')
        gc = ls[1].split('|')[0]
        gc_is_sets[gc].add(ls[1])

print('DATASET_SIMPLEBAR')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tIS Counts')
print('DATA')
for gc in gc_is_sets:
    if gc in mapping:
        print(mapping[gc] + '\t' + str(len(gc_is_sets[gc])))
