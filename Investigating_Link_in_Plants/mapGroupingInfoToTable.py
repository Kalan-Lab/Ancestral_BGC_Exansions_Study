import os
import sys

early_diverging_clades = ['Marchantiophyta', 'Bryophyta', 'Charophyceae', 'Klebsormidiales']

early_diverging_genera = set([])
with open('Lineage_Info.txt') as olif:
    for line in olif:
        line = line.strip()
        ls = line.split('\t')
        for edc in early_diverging_clades:
            if edc in ls[1]:
                early_diverging_genera.add(ls[1].split(';')[-1])

with open('Data_from_PlantiSMASH_Paper.txt') as odf:
    for i, line in enumerate(odf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: 
            print(line + '\tgroup')
            continue
        genus = ls[2].split()[0]
        group = 'Vascular Land Plants'
        if genus in early_diverging_genera:
            group = 'Early Diverging Streptophyta'
        print('\t'.join(ls + [group]))

