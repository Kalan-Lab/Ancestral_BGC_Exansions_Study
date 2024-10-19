import os
import sys
from collections import defaultdict
import statistics

clade_map = {}
clade_dir = 'Clades/'

for cf in os.listdir(clade_dir):
    clade_file = clade_dir + cf
    clade = cf.split('.txt')[0]
    with open(clade_file) as ocf:
        for line in ocf:
            line = line.strip()
            clade_map[line] = clade

# gca     name    bgc     bgc_length      completeness    putative        metallophore-like       nrps-like       pks-like        terpene-like    product-types   bgc_gbk
cats = ['Metallophore', 'NRPS', 'PKS', 'Terpene',  'Other']
clade_type_lengths = defaultdict(lambda: defaultdict(list))
with open('BGC_Details_Fungal.txt') as obdf:
    for i, line in enumerate(obdf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        clade = None
        if ls[1] in clade_map:
            clade = clade_map[ls[1]]
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
        if len(clade_type_lengths[clade][cat]) >= 1:
            clade_type_length_medians[clade][cat] = statistics.median(clade_type_lengths[clade][cat])
        else:
            clade_type_length_medians[clade][cat] = 1e-100

tot_clade_counts = defaultdict(lambda: defaultdict(float))
full_clade_type_counts = defaultdict(lambda: defaultdict(float))
strict_clade_type_counts = defaultdict(lambda: defaultdict(float))
with open('BGC_Details_Fungal.txt') as obdf:
    for i, line in enumerate(obdf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        clade = None
        if ls[1] in clade_map:
            clade = clade_map[ls[1]]
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

print('group\tcat\t' + '\t'.join(cats))
for clade in full_clade_type_counts:
    row1 = [clade, 'full']
    row2 = [clade, 'strict'] 
    for cat in cats:
        row1.append(str(100.0*(full_clade_type_counts[clade][cat]/tot_clade_counts[clade]['full'])))
        row2.append(str(100.0*(strict_clade_type_counts[clade][cat]/tot_clade_counts[clade]['strict'])))
    print('\t'.join(row1))
    print('\t'.join(row2))
