import os
import sys
from ete3 import Tree


fields = ['MIBiG', 'MIBiG (benzenic)', 'MIBiG (naphthalenic)', 'Mycobacteriales', 'Streptomycetales', 'Streptosporangiales']
field_colors = ['#4287f5', '#97b6e6', '#122e59', '#FF66CC', '#9966FF', '#9966FF']

print('DATASET_PIECHART')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tMIBiG')
print('FIELD_LABELS\t' + '\t'.join(fields))
print('FIELD_COLORS\t' + '\t'.join(field_colors))
print('DATA')

taxa_file = 'GCA_to_Taxonomy.txt'
mibig_class_file = 'MIBiG_Classification.txt'

mibig = {}
with open(mibig_class_file) as omcf:
    for line in omcf:
        line = line.strip('\n')
        ls = line.split('\t')
        mibig[ls[0]] = ls[1]

taxa = {}
with open(taxa_file) as otf:
    for line in otf:
        line = line.strip()
        ls = line.split('\t')
        if not ';o__' in line: continue
        order = ls[1].split(';o__')[1].split(';')[0]
        if order in fields:
            taxa[ls[0]] = order

tre = 'RifK_Homologs_including_MIBiG_and_Bacillota.tre'
t = Tree(tre)

for n in t.traverse('postorder'):
    if n.is_leaf():
        name = n.name
        row = []
        for f in fields:
            if 'MIBiG' in name:
                mibig_group = mibig[name]
                if mibig_group == f:
                    row.append('1')
                else:
                    row.append('0')
            else:
                gca = name.split('|')[0]
                if gca in taxa and taxa[gca] == f:
                    row.append('1')
                else:
                    row.append('0')
        if sum([int(x) for x in row]) > 0 and 'MIBiG' in name: 
            print('\t'.join([name, '1', '1'] + row))
