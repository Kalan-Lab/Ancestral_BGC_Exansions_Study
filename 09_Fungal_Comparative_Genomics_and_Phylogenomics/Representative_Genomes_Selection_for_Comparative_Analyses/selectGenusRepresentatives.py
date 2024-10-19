import os
import sys
from collections import defaultdict
from operator import itemgetter

as_file = 'AntiSMASH_Stats.txt'

genus_gcas = defaultdict(list)
with open(as_file) as oaf:
    for i, line in enumerate(oaf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        genus = ls[1].split('_')[0]
        genus_gcas[genus].append([ls[0], float(ls[4])])

for genus in genus_gcas:
    for i, go in enumerate(sorted(genus_gcas[genus], key=itemgetter(1), reverse=True)):
        if i == 0:
            print(go[0] + '\t' + str(go[1]))
