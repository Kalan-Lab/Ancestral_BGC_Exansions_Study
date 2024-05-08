import os
import sys

dikarya = set([])
with open('Dikarya_Names.txt') as odnf:
    for line in odnf:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        dikarya.add(gca)

with open('AntiSMASH_Stats.Updated_WithoutContaminants.txt') as oaf:
    for i, line in enumerate(oaf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        if not ls[0] in dikarya: continue
        bos = float(ls[-1])/1000000.0
        print(ls[0] + '\t' + str(bos))
