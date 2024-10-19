import os
import sys

as_tsv = 'AntiSMASH_Stats.txt'

with open(as_tsv) as oat:
    for i, line in enumerate(oat):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        gca = ls[0]
        bgc_sum = float(ls[4])/1000000
        print(gca  + '\t' + str(bgc_sum))
