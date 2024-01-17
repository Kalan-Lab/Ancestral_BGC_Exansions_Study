import os
import sys
from collections import defaultdict
from ete3 import Tree
from scipy import stats

clade_1_file = 'Clade-1.txt'
clade_2_file = 'Clade-2.txt'
bgc_file = 'BGC_Gene_OGs.txt'

clade_1 = set([])
clade_2 = set([])

with open(clade_1_file) as oc1f:
    for line in oc1f:
        line = line.strip()
        clade_1.add(line)

with open(clade_2_file) as oc2f:
    for line in oc2f:
        line = line.strip()
        clade_2.add(line)

og_clade_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
og_genes_clade = defaultdict(dict)
with open(bgc_file) as obf:
    for line in obf:
        line = line.strip()
        ls = line.split('\t')
        og = ls[3]
        if og == 'NA': continue
        lt = ls[2]
        genus = ls[1]
        if genus in clade_1:
            og_clade_counts[og]['clade_1'][genus] += 1
            og_genes_clade[og][lt] = 'clade_1'
        elif genus in clade_2:
            og_clade_counts[og]['clade_2'][genus] += 1
            og_genes_clade[og][lt] = 'clade_2'
        else:
            og_clade_counts[og]['other'][genus] += 1
            og_genes_clade[og][lt] = 'other'

for og in og_clade_counts:
    ogc1 = sum(og_clade_counts[og]['clade_1'].values())
    ogc2 = sum(og_clade_counts[og]['clade_2'].values())
    tot = ogc1 + ogc2
    if tot >= 2 and ogc1 >= 1 and ogc2 >= 1:
        print(og + '\t' + str(ogc1) + '\t' + str(ogc2) + '\t' + ', '.join([x[0] for x in og_clade_counts[og]['clade_1'].items()]) + '\t' + ', '.join([x[0] for x in og_clade_counts[og]['clade_2'].items()]) + '\t' + ', '.join([x[0] for x in og_clade_counts[og]['other'].items()]))
