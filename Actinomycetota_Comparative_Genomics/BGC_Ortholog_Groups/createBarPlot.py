import os
import sys
from collections import defaultdict

heatmap_track_file = 'BGC_OG_Heatmap.iTol.txt'
class_file = 'Actinomycetia.txt'
og_to_prod_file = 'OG_Info_for_SupTable.txt'

og_product_category = {}
with open(og_to_prod_file) as oopf:
    for line in oopf:
        line = line.strip() 
        ls = line.split('\t')
        og = ls[0]
        prods = ls[-2]
        prod = []
        max_count = 0.0
        sum_count = 0.0
        for p in prods.split('; '):
            prod.append(p.split('=')[0])
            sum_count += float(p.split('=')[1])
            if max_count < float(p.split('=')[1]):
                max_count = float(p.split('=')[1])
        
        max_freq = max_count/sum_count
        product = 'NA'
        product_string = ' '.join(prod) + ' '
        if len(prod) == 1 or max_freq >= 0.80:
            product = '(largely) specific BGC type'
        elif len(prod) == 2 and 'NRPS-like' in prod and 'NRPS' in prod:
            product = '(largely) specific BGC type'
        elif len(prod) >= 2:
            if 'NRPS ' in product_string and 'PKS' in product_string:
                product = 'multi-type (with NRPS & PKS)'
            elif 'NRPS ' in product_string:
                product = 'multi-type (with NRPS)'
            elif 'PKS' in product_string:
                product = 'multi-type (with PKS)'
            else:
                product = 'multi-type'

        og_product_category[og] = product

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

conserved_bgc = set([])
conserved_other = set([])
completeley_new = set([])
for og in ogs:
    #other_total = bgc_context[og]['Other'] + other_context[og]['Other']
    #actino_total = bgc_context[og]['Actinomycetia'] + other_context[og]['Actinomycetia']
    if bgc_context[og]['Other'] > 0:
        conserved_bgc.add(og)
    elif other_context[og]['Other'] > 0:
        conserved_other.add(og)
    else:
        completeley_new.add(og)

#print('General\t' + str(len(conserved)))
#print('Actinomycetia Class-Specific\t' + str(len(class_specific)))

comp_new_counts = defaultdict(int)
cons_bgc_counts = defaultdict(int)
cons_oth_counts = defaultdict(int)
for og in completeley_new:
    product = og_product_category[og]
    comp_new_counts[product] += 1

for og in conserved_bgc:
    product = og_product_category[og]
    cons_bgc_counts[product] += 1

for og in conserved_other:
    product = og_product_category[og]
    cons_oth_counts[product] += 1

for annot in comp_new_counts:
    print('Actinomycetia specific\t' + annot + '\t' + str(comp_new_counts[annot]))

for annot in cons_bgc_counts:
    print('Ancestral - in a BGC context\t' + annot + '\t' + str(cons_bgc_counts[annot]))

for annot in cons_oth_counts:
    print('Ancestral - not in a BGC context\t' + annot + '\t' + str(cons_oth_counts[annot]))
