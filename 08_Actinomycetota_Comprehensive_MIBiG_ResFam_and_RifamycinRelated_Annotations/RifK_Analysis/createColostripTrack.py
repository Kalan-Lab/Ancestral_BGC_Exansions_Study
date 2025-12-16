import os
import sys
from ete3 import Tree


print('DATASET_COLORSTRIP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tTaxa')
print('DATA')

taxa_file = 'GCA_to_Taxonomy.txt'

taxa = {}
with open(taxa_file) as otf:
    for line in otf:
        line = line.strip()
        ls = line.split('\t')
        if not ';o__' in line: continue
        order = ls[1].split(';o__')[1].split(';')[0]
        taxa[ls[0]] = order

tre = 'RifK_Homologs_including_MIBiG_and_Bacillota.tre'
t = Tree(tre)

for n in t.traverse('postorder'):
    if n.is_leaf():
        name = n.name
        gca = name.split('|')[0]
        if gca in taxa:
            tax_col = None
            if taxa[gca] == 'Mycobacteriales':
                tax_col = '#FF66CC'
            elif taxa[gca] == 'Streptomycetales' or taxa[gca] == 'Streptosporangiales':
                tax_col = '#9966FF'
            if tax_col != None:
                print(name + '\t' + tax_col + '\t' + taxa[gca])
