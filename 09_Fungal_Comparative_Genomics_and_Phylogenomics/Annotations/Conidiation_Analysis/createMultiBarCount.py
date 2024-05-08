import os
import sys
from collections import defaultdict

gca_to_id = {}
with open('../AntiSMASH_Stats.Updated_WithoutContaminants.txt') as oaf:
    for i, line in enumerate(oaf):
        line = line.strip()
        ls = line.split('\t')
        gca_to_id[ls[0]] = ls[1]

select_ogs = set([])
select_hogs = set([])
with open('OG_Mappings.txt') as ogmf:
    for line in ogmf:
        line = line.strip()
        ls = line.split('\t')
        select_ogs.add(ls[1])
        select_hogs.add(ls[2])

gca_og_counts = defaultdict(int)
gca_hog_counts = defaultdict(int)

gcas = []
with open('Orthogroups.tsv') as ogf:
    for i, line in enumerate(ogf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            gcas = ls[1:]
        else:
            og = ls[0]
            if not og in select_ogs: continue
            for j, gs in enumerate(ls[1:]):
                gca = gcas[j]
                if gs.strip() != '':
                    gca_og_counts[gca] += 1 # len(gs.split(', '))

gcas = []
with open('N0.tsv') as onf:
    for i, line in enumerate(onf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: 
            gcas = ls[3:]
        else:
            hog = ls[0]
            if not hog in select_hogs: continue
            for j, gs in enumerate(ls[3:]):
                gca = gcas[j]
                if gs.strip() != '':
                    gca_hog_counts[gca] += 1# len(gs.split(', '))


print('DATASET_SIMPLEBAR')
#print('DATASET_MULTIBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tConidiation-associated homologs')
print('COLOR\t#000000')
#print('FIELD_LABELS\tHOGs\tNon-HOG_OGs')
#print('FIELD_COLORS\t#8a0707\t#e65e5e')
print('DATA')
for gca in gcas:
    hog_count = gca_hog_counts[gca]
    non_hog_count = gca_og_counts[gca] - hog_count
    print(gca_to_id[gca] + '\t' + str(hog_count))# + '\t' + str(non_hog_count))
