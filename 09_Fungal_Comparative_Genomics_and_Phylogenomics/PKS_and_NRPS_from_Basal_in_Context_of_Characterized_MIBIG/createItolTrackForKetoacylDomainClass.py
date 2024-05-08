import os
import sys
from collections import defaultdict
from ete3 import Tree

domain_class = {'type I iterative cis-AT': '#f5c462', 'type I modular cis-AT': '#c4932f', 'type I trans-AT': '#ab7509', 'other/unclassified': '#818285'}

seq_to_dc = defaultdict(lambda: 'other/unclassified')
with open('NAPDOS2_KS_Data.txt') as of:
    for i, line in enumerate(of):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        if ls[6] in domain_class.keys():
            seq_to_dc['_'.join(ls[0].split('_')[:-2])] = ls[6] 

print('DATASET_COLORSTRIP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tDomType')
print('DATA')

t = Tree('KetoacylSynthase_Domain_Sequences.tre')
for n in t.traverse('postorder'):
    if n.is_leaf():
        nn = n.name
        dc = seq_to_dc[nn]
        dc_color = domain_class[dc]
        print(nn + '\t' + dc_color + '\t' + dc)
