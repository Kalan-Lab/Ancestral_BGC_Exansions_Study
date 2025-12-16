import os
import sys

ga_cutoff = 22.50

gather_set = set([])
with open('../HI_Full_Hits.Domain.Cut_GA.txt') as ohif:
    for line in ohif:
        line = line.strip()
        if line.startswith('#'): continue
        ls = line.split()
        name = ls[0] + '|' + ls[17] + '|' + ls[18]
        gather_set.add(name)

print('domain_evalue\tdomain_score\tmeets_ga')
with open('HI_Full_Hits.Domain.Evalue_10.txt') as ohif:
    for line in ohif:
        line = line.strip()
        if line.startswith('#'): continue
        ls = line.split()
        if ls[3] != 'HET': continue
        name = ls[0] + '|' + ls[17] + '|' + ls[18]
        hit_flag = "False"
        if name in gather_set:
           hit_flag = "True"
        print(name + '\t' + ls[12] + '\t' + ls[13] + '\t' + hit_flag)
