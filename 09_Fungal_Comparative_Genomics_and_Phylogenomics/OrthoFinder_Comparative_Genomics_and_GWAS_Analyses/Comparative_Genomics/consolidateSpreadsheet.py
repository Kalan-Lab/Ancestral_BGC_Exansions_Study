import os
import sys
from collections import defaultdict

eggnog_annots = defaultdict(lambda: [0.0, (['NA']*20)])
header = []
with open('out.emapper.annotations') as oeaf:
    for line in oeaf:
        line = line.strip('\n')
        ls = line.split('\t')
        if not line.startswith('#'):
            score = float(ls[3])
            og = ls[0].split('|')[0]
            if score > eggnog_annots[og][0]:
                eggnog_annots[og] = [score, ls]
        elif line.startswith('#query'):
            header = [ls[0][1:]] + ls[1:]

with open('Comparative_Results.txt') as oshf:
    for i, line in enumerate(oshf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            print('\t'.join(ls + header))
        else:
            og = ls[0]
            print('\t'.join(ls + eggnog_annots[og][1]))
