import os
import sys

gca_to_name = {}
for line in open('../FungiWide_Genome_GCAS.txt'):
    line = line.strip()
    ls = line.split('\t')
    gca_to_name[ls[-1].split('.')[0]] = ls[0]

for i, line in enumerate(open('Concatenated_Alignment.phy')):
    line = line.strip()
    if i == 0: 
        print(line)
    else:
        aid = line.split()[0]
        gca = aid.split('.')[0]
        name = gca_to_name[gca]
        print(line.replace(aid, name))
        
