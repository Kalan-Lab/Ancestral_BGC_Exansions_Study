import os
import sys
from collections import defaultdict

valid_architectures = set(['Het-C', 'ATP-cone|Ribonuc_red_lgN|Ribonuc_red_lgC', 'HMG_box', 'MAT1-1-2', 'Patatin|NB-ARC', 'HeLo|HET-s_218-289', 'HET|NACHT|WD40', 'SET|Rubis-subs-bind|zf-MYND', 'Peptidase_S8', 'GLTP', 'RVT_2', 'SERF-like_N', 'CYSTM', 'Ccdc124', 'Fructosamin_kin', 'Patatin|TPR_10'])

gca_names = {}
with open('../../Overview_File_with_Corrected_N50s.txt') as ovf:
    for line in ovf:
        line = line.strip()
        ls = line.split('\t')
        gca_names[ls[0]] = ls[1]

gca_het_hits = defaultdict(set)
gca_other_hits = defaultdict(set)
gcas = set([])
with open('Domain_Architectures_of_Proteins_with_OneOrMoreHIAssociatedDoms.txt') as ohh:
    for line in ohh:
        line = line.strip()
        ls = line.split('\t')
        gene = ls[0]
        gca = ls[0].split('|')[0]
        if ls[1] == 'HET':
            gca_het_hits[gca_names[gca]].add(gene)
        elif ls[1] in valid_architectures:
            gca_other_hits[gca_names[gca]].add(gene)
        gcas.add(gca_names[gca])

print('DATASET_MULTIBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tHet/vegetative incompatability genes')
print('FIELD_COLORS\t#7a3e56\t#aa5ed6')
print('FIELD_LABELS\tSingleton_HET\tOther')
print('DATA')

for gca in gcas:
    count_with_het = len(gca_het_hits[gca])
    count_remaining = len(gca_other_hits[gca])
    print(gca + '\t' + str(count_with_het) + '\t' + str(count_remaining))
