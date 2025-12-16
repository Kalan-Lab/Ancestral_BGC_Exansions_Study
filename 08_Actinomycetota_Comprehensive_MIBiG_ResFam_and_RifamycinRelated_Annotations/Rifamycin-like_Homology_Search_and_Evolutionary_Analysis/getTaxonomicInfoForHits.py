import os
import sys

gca_class = {}
gca_fam = {}
with open('../GCA_to_Taxonomy.txt') as ogt:
    for line in ogt:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('.')[0]
        taxclass = ls[1].split('c__')[1].split(';o__')[0]
        taxfam = ls[1].split('f__')[1].split(';g__')[0]
        gca_class[gca] = taxclass
        gca_fam[gca] = taxfam

clade1 = set([])
clade2 = set([])
with open('Clade_1_Families.txt') as oc1f:
    for line in oc1f:
        line = line.strip()
        clade1.add(line)
with open('Clade_2_Families.txt') as oc2f:
    for line in oc2f:
        line = line.strip()
        clade2.add(line)

with open('Tiny_AAI_Plot_Data.txt') as otapd:
    for i, line in enumerate(otapd):
        line = line.strip()
        if i == 0:
            print('ClassGroup\t' + line)
        else:
            gca = line.split('.')[0]
            gclass = gca_class[gca]
            gfam = gca_fam[gca]
            fam_group = 'Other Actinomycetota'
            if gclass == 'Actinomycetia':
                fam_group = 'Actinomycetia'
                if gfam in clade1:
                    fam_group = 'Actinomycetia - Clade-1'
                elif gfam in clade2:
                    fam_group = 'Actinomycetia - Clade-2'
            print(fam_group + '\t' + line)
