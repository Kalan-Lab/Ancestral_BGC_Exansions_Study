import os
import sys
from collections import defaultdict
from scipy import stats
import numpy 

heatmap_track_file = 'BGC_OG_Heatmap.iTol.txt'
class_file = 'Actinomycetia.txt'
og_to_prod_file = 'BGC_Gene_OGs.With_Product_Info.txt'

products = defaultdict(lambda: defaultdict(float))
with open(og_to_prod_file) as of:
    for line in of:
        line = line.strip()
        ls = line.split('\t')
        og = ls[3]
        prod = ls[-1]
        for p in prod.split(' '):
            products[og][p] += 1.0/float(len(prod.split(' ')))

actinos = set([])
with open(class_file) as ocf:
    for line in ocf:
        line = line.strip()
        actinos.add(line)

bgc_context = defaultdict(lambda: defaultdict(int))
other_context = defaultdict(lambda: defaultdict(int))
ogs = []
with open(heatmap_track_file) as of:
    for i, line in enumerate(of):
        line = line.strip()
        ls = line.split('\t')
        if i == 4:
            ogs = ls[1:]
        elif i >= 6:
            leaf = ls[0]
            group = 'Other' 
            if leaf in actinos:
                group = 'Actinomycetia'
            for j, val in enumerate(ls[1:]):
                og = ogs[j]
                if val == '1':
                    other_context[og][group] += 1
                elif val == '2':
                    bgc_context[og][group] += 1

for og in ogs:
    prods = products[og]
    bgc_type_string_list = []
    for p in prods:
        bgc_type_string_list.append(p + '=' + str(prods[p]))
    bgc_type_string = '; '.join(bgc_type_string_list)

    bgc_actino = bgc_context[og]['Actinomycetia'] 
    other_actino = other_context[og]['Actinomycetia']
    bgc_other = bgc_context[og]['Other']
    other_other = other_context[og]['Other']
    stat, pval = stats.fisher_exact([[bgc_actino, bgc_other], [other_actino, other_other]])
    count_info = '(' + str(bgc_actino) + ', ' + str(other_actino) + ', ' + str(bgc_other) + ', ' + str(other_other) + ')'
    group = 'Actinomycetia specific'
    if bgc_other > 0:
        group = 'Ancestrally conserved; found in a BGC context'
    elif other_other > 0:
        group = 'Ancestrally conserved; not found in a BGC context'

    print(og + '\t' + group + '\t' + bgc_type_string + '\t' + str(pval) + ' ' + count_info)

