import os
import sys
from collections import defaultdict

heatmap_track_file = 'BGC_OG_Heatmap.iTol.txt'

print('DATASET_MULTIBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tTotals')
print('COLOR\t#000000')
print('FIELD_LABELS\t' + '\t'.join(['BGC Context', 'Other Context']))
print('FIELD_COLORS\t' + '\t'.join(['#000000', '#AFABAB']))
print('DATA')

bgc_context = defaultdict(int)
other_context = defaultdict(int)
all_leafs = set([])
with open(heatmap_track_file) as of:
    for i, line in enumerate(of):
        if i >= 6:
            line = line.strip()
            ls = line.split('\t')
            leaf = ls[0]
            for val in ls[1:]:
                if val == '1':
                    other_context[leaf] += 1
                elif val == '2':
                    bgc_context[leaf] += 1
            all_leafs.add(leaf)

for leaf in all_leafs:
    print(leaf + '\t' + str(bgc_context[leaf]) + '\t' + str(other_context[leaf]))    

