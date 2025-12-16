import os
import sys
from collections import defaultdict

gca_names = {}
with open('../../Overview_File.txt') as ovf:
    for line in ovf:
        line = line.strip()
        ls = line.split('\t')
        gca_names[ls[0]] = ls[1]

gca_het_hits = defaultdict(set)
gca_nac_hits = defaultdict(set)
gcas = set([])
with open('HI_Hits.Domain.Cut_TC.txt') as ohh:
    for line in ohh:
        if line.startswith('#'): continue
        line = line.strip()
        ls = line.split()
        if ls[3] == 'HET':
            gene = ls[0]
            gca = gene.split('|')[0]
            gca_het_hits[gca_names[gca]].add(gene)
            gcas.add(gca_names[gca])
        elif ls[3] == 'NACHT':
            gene = ls[0]
            gca = gene.split('|')[0]
            gca_nac_hits[gca_names[gca]].add(gene)
            gcas.add(gca_names[gca])

print('DATASET_MULTIBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tHet/vegetative incompatability genes')
print('FIELD_COLORS\t#7a3e56\t#aa5ed6')
print('FIELD_LABELS\tHET\tNACHT_Only')
print('DATA')

for gca in gcas:
    count_with_het = len(gca_het_hits[gca])
    count_remaining = len(gca_nac_hits[gca].difference(gca_het_hits[gca]))
    print(gca + '\t' + str(count_with_het) + '\t' + str(count_remaining))
