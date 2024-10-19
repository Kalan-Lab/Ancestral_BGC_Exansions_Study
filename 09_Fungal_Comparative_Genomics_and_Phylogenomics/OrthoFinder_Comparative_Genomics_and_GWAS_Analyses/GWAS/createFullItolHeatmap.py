import os
import sys
from collections import defaultdict

result_file = 'GWAS_Results.txt'
og_file = 'Orthogroups.tsv'
ov_file = '../../Overview_File.txt'

naming = {}
with open(ov_file) as ovf:
    for i, line in enumerate(ovf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        naming[ls[0]] = ls[1]

select = set([])
og_order = []
og_named_order = []
with open(result_file) as osf:
    for line in osf:
        og = line.strip().split('\t')[0]
        select.add(og)
        og_order.append(og)
        og_named_order.append(og)

genomes = []
gca_og_carriage = defaultdict(set)
with open(og_file) as oof:
    for i, line in enumerate(oof):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            genomes = ls[1:]
        else:
            og = ls[0]
            if og in select:
                for j, gs in enumerate(ls[1:]):
                    gca = genomes[j]
                    if gs.strip() != '':
                        gca_og_carriage[gca].add(og)
                    
print('DATASET_HEATMAP')
print('SEPARATOR TAB')
print('DATASET_LABEL\tGWAS Results')
print('COLOR\t#000000')
print('FIELD_LABELS\t' + '\t'.join(og_named_order))
print('DATA')
for gca in gca_og_carriage:
    printlist = [naming[gca]]
    for og in og_order:
        if og in gca_og_carriage[gca]:
            printlist.append('1')
        else:
            printlist.append('0')
    print('\t'.join(printlist))
