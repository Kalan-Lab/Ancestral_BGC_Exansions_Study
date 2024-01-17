import os
import sys
from collections import defaultdict
from ete3 import Tree

t = Tree('Homologs_to_Sv_HelR_Full.msa.trimal_strict.tre')
family_info_file = 'Family_to_Type_Mapping.txt'

family_to_info = {}
with open(family_info_file) as ofif:
    for line in ofif:
        line = line.strip()
        ls = line.split('\t')
        family_to_info[ls[0]] = [ls[1], ls[2]]


print('DATASET_COLORSTRIP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tFamily_Grouping')
print('DATA')

for n in t.traverse('postorder'):
    if n.is_leaf:
        leafname = n.name
        family = leafname.split('|')[0]
        if family in family_to_info:
            print(leafname + '\t' + family_to_info[family][1] + '\t' + family_to_info[family][0])
