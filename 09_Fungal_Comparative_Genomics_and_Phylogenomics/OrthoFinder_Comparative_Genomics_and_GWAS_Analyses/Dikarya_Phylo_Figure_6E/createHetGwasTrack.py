import os
import sys
from ete3 import Tree
from collections import defaultdict

het_hogs = set([])
with open('HET_GWAS_HOGs.txt') as ohgh:
    for i, line in enumerate(ohgh):
        line = line.strip()
        het_hogs.add(line) 

gcas = []
gca_hogs = defaultdict(set)
with open('../GWAS_Analysis/Orthogroup_Matrix.txt') as ogom:
    for i, line in enumerate(ogom):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            gcas = ls[1:]
            continue
        else:
            hog = ls[0].split('N0.')[1]
            if hog in het_hogs:
                for j, val in enumerate(ls[1:]):
                    if val == '1':
                        gca = gcas[j]
                        gca_hogs[gca].add(hog)

t = Tree('Fungi_Wide_Tree_with_Outgroups_pared.tre') 

print('DATASET_HEATMAP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tGWAS HET HOGs')
print('FIELD_LABELS\t' + '\t'.join(sorted(het_hogs)))
print('DATA')

for n in t.traverse('postorder'):
    if n.is_leaf():
        gca = '_'.join(n.name.split('_')[-2:])
        printlist = [n.name]
        for hh in sorted(het_hogs):
            if hh in gca_hogs[gca]:
                printlist.append('1')
            else:
                printlist.append('0')
        print('\t'.join(printlist))
