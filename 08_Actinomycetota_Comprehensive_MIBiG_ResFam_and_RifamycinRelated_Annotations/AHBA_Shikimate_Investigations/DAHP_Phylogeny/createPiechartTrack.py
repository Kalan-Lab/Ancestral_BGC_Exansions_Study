import os
import sys
from ete3 import Tree

print('DATASET_PIECHART')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tRifH_Homologs')
print('FIELD_LABELS\tAN')
print('FIELD_COLORS\t#CC66E6')
print('DATA')

t = Tree('DAHP_synthase_with_RifH_Homologs.tre')

for n in t.traverse('postorder'):
    if n.is_leaf():
        nname = n.name
        if 'MIBiG' in nname:
            print('\t'.join([nname, '1', '1', '1']))
