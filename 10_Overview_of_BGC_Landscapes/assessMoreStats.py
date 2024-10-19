import os
import sys
from collections import defaultdict
import statistics

bacterial_orders = {'Cyanobacteriales': 'Cyanobacteriota', 'Mycobacteriales': 'Actinomycetota', 'Streptosporangiales': 'Actinomycetota', 'Streptomycetales': 'Actinomycetota', 'Myxococcales': 'Myxococcota'}

pezi = set([])
agari = set([])
with open('BGC_Enriched_Pezizomycotina.txt') as obepf:
    for line in obepf:
        line = line.strip()
        pezi.add(line)

with open('Agaricomycetes.txt') as oaf:
    for line in oaf:
        line = line.strip()
        agari.add(line)

# gca     name    bgc     bgc_length      completeness    putative        metallophore-like       nrps-like       pks-like        terpene-like    product-types   bgc_gbk
cats = ['Metallophore', 'NRPS', 'PKS', 'Terpene',  'Other']
clade_type_lengths = defaultdict(lambda: defaultdict(list))
with open('BGC_Details_Fungal.txt') as obdf:
    for i, line in enumerate(obdf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        clade = None
        if ls[1] in pezi:
            clade = 'BGC_Enriched_Pezizomycotina'
        elif ls[1] in agari:
            clade = 'Agaricomycetes'
        else: continue
        bgc_length = float(ls[3])
        completeness = ls[4]
        if completeness == 'True':
            none_true = True
            for j, boolean in enumerate(ls[6:-2]):
                bgc_cat = cats[j]
                if boolean == 'True':
                    none_true = False
                    clade_type_lengths[clade][bgc_cat].append(bgc_length)
            if none_true:
                bgc_cat = 'Other'
                clade_type_lengths[clade][bgc_cat].append(bgc_length)

with open('BGC_Details_Bacterial.txt') as obdf:
    for i, line in enumerate(obdf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        clade = None
        order = ls[1].split('|')[0]
        if order in bacterial_orders:
            clade = bacterial_orders[order]
        else: continue
        bgc_length = float(ls[3])
        completeness = ls[4]
        if completeness == 'True':
            none_true = True
            for j, boolean in enumerate(ls[6:-2]):
                bgc_cat = cats[j]
                if boolean == 'True':
                    none_true = False
                    clade_type_lengths[clade][bgc_cat].append(bgc_length)
            if none_true:
                bgc_cat = 'Other'
                clade_type_lengths[clade][bgc_cat].append(bgc_length)

clade_type_length_medians = defaultdict(lambda: defaultdict(float))
for clade in clade_type_lengths:
    for cat in cats:
        clade_type_length_medians[clade][cat] = statistics.median(clade_type_lengths[clade][cat])

tot_clade_counts = defaultdict(lambda: defaultdict(float))
full_clade_type_counts = defaultdict(lambda: defaultdict(float))
strict_clade_type_counts = defaultdict(lambda: defaultdict(float))
multi_type_counts = defaultdict(int)
tot_bgc_counts = defaultdict(int)
with open('BGC_Details_Fungal.txt') as obdf:
    for i, line in enumerate(obdf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        clade = None
        if ls[1] in pezi:
            clade = 'BGC_Enriched_Pezizomycotina'
        elif ls[1] in agari:
            clade = 'Agaricomycetes'
        else: continue
        bgc_length = float(ls[3])
        completeness = ls[4]
        putative = ls[5]
        none_true = True
        for j, boolean in enumerate(ls[6:-2]):
            bgc_cat = cats[j]
            if boolean == 'True':
                none_true = False
                if completeness == 'True':
                    full_clade_type_counts[clade][bgc_cat] += 1
                    if putative == 'False':
                        strict_clade_type_counts[clade][bgc_cat] += 1
                else:
                    count_prop = min([1.0, bgc_length/clade_type_length_medians[clade][cat]])
                    full_clade_type_counts[clade][bgc_cat] += count_prop
                    if putative == 'False':
                        strict_clade_type_counts[clade][bgc_cat] += count_prop
        if none_true:
            bgc_cat = 'Other'
            if completeness == 'True':
                full_clade_type_counts[clade][bgc_cat] += 1
                if putative == 'False':
                    strict_clade_type_counts[clade][bgc_cat] += 1
            else:
                count_prop = min([1.0, bgc_length/clade_type_length_medians[clade][cat]])
                full_clade_type_counts[clade][bgc_cat] += count_prop
                if putative == 'False':
                    strict_clade_type_counts[clade][bgc_cat] += count_prop
        if completeness == 'True':
            tot_clade_counts[clade]['full'] += 1
            if putative == 'False':
                tot_clade_counts[clade]['strict'] += 1
        else:
            count_prop = min([1.0, bgc_length/clade_type_length_medians[clade][cat]])
            tot_clade_counts[clade]['full'] += count_prop
            if putative == 'False':
                tot_clade_counts[clade]['strict'] += count_prop

        if completeness =='True':
            tot_bgc_counts[clade] += 1
            products = set(ls[-2].split('|'))
            if len(products) > 1 and (not (len(products) == 2 and 'NRPS-like' in products and 'NRPS' in products)):
                multi_type_counts[clade] += 1

with open('BGC_Details_Bacterial.txt') as obdf:
    for i, line in enumerate(obdf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        clade = None
        order = ls[1].split('|')[0]
        if order in bacterial_orders:
            clade = bacterial_orders[order]
        else: continue
        bgc_length = float(ls[3])
        completeness = ls[4]
        putative = ls[5]
        none_true = True
        for j, boolean in enumerate(ls[6:-2]):
            bgc_cat = cats[j]
            if boolean == 'True':
                none_true = False
                if completeness == 'True':
                    full_clade_type_counts[clade][bgc_cat] += 1
                    if putative == 'False':
                        strict_clade_type_counts[clade][bgc_cat] += 1
                else:
                    count_prop = min([1.0, bgc_length/clade_type_length_medians[clade][cat]])
                    full_clade_type_counts[clade][bgc_cat] += count_prop
                    if putative == 'False':
                        strict_clade_type_counts[clade][bgc_cat] += count_prop
        if none_true:
            bgc_cat = 'Other'
            if completeness == 'True':
                full_clade_type_counts[clade][bgc_cat] += 1
                if putative == 'False':
                    strict_clade_type_counts[clade][bgc_cat] += 1
            else:
                count_prop = min([1.0, bgc_length/clade_type_length_medians[clade][cat]])
                full_clade_type_counts[clade][bgc_cat] += count_prop
                if putative == 'False':
                    strict_clade_type_counts[clade][bgc_cat] += count_prop
        if completeness == 'True':
            tot_clade_counts[clade]['full'] += 1
            if putative == 'False':
                tot_clade_counts[clade]['strict'] += 1
        else:
            count_prop = min([1.0, bgc_length/clade_type_length_medians[clade][cat]])
            tot_clade_counts[clade]['full'] += count_prop
            if putative == 'False':
                tot_clade_counts[clade]['strict'] += count_prop
        
        if completeness == 'True':
            tot_bgc_counts[clade] += 1
            products = set(ls[-2].split('|'))
            if len(products) > 1 and (not (len(products) == 2 and 'NRPS-like' in products and 'NRPS' in products)):
                multi_type_counts[clade] += 1

"""
print('group\tbgc_category\tfull_count\tstrict_count\tmedian_length')
for clade in full_clade_type_counts:
    for cat in cats:
        print('\t'.join([str(x) for x in [clade, cat, full_clade_type_counts[clade][cat], strict_clade_type_counts[clade][cat], clade_type_length_medians[clade][cat]]]))
"""

print('clade\ttot_bgc_count\tmulti_type_proportion')
for clade in tot_bgc_counts:
    print(clade + '\t' + str(tot_bgc_counts[clade]) + '\t' + str(100.0*(multi_type_counts[clade]/tot_bgc_counts[clade])))
