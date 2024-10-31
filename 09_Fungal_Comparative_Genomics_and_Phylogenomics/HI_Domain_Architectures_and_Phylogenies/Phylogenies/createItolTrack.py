import os
import sys
from ete3 import Tree

agar = set([])
pezi = set([])

with open('../Clades/Agaricomycetes.txt') as oaf:
    for line in oaf:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:]).split('.')[0].replace('GCF_', 'GCA_')
        agar.add(gca)

with open('../Clades/Pezizomycotina.txt') as oaf:
    for line in oaf:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:]).split('.')[0].replace('GCF_', 'GCA_')
        pezi.add(gca)

print('DATASET_COLORSTRIP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tClades')
print('DATA')
t = Tree(sys.argv[1])

for n in t.traverse('postorder'):
    if n.is_leaf():
        nname = n.name.split('|')[0]
        clade = 'other'
        color = '#8c8a89'
        if nname in agar:
            clade = 'Agaricomycetes'
            color = '#FF5D5D'
        elif nname in pezi:
            clade = 'Pezizomycotina'
            color = '#9ECDF8'
        print(n.name + '\t' + color + '\t' + clade) 
