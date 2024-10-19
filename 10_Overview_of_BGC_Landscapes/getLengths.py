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
print('clade\tbgc_cat\tbgc_length')
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
                    print(clade + '\t' + bgc_cat + '\t' + str(bgc_length))
            if none_true:
                bgc_cat = 'Other'
                print(clade + '\t' + bgc_cat + '\t' + str(bgc_length))

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
                    print(clade + '\t' + bgc_cat + '\t' + str(bgc_length))
            if none_true:
                bgc_cat = 'Other'
                print(clade + '\t' + bgc_cat + '\t' + str(bgc_length))
