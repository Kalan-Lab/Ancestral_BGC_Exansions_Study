import os
import sys
from collections import defaultdict

og_hits = defaultdict(set)
og_total = defaultdict(set)

hit_cds = set([])
with open('Heterokaryon_Search_All_Prots.txt') as of:
    for line in of:
        line = line.strip()
        if line.startswith('#'): continue
        ls = line.split()
        evalue = float(ls[4])
        if evalue < 1e-3:
            hit_cds.add(ls[0])
            #gca_hits[ls[0].split('|')[0]] += 1

gcas = []
with open(sys.argv[1]) as of:
    for i, line in enumerate(of):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            gcas = ls[1:]
            continue
        og = ls[0]
        for j, cdss in enumerate(ls[1:]):
            gca = gcas[j]
            for cds in cdss.split(', '):
                if cds.strip() != '' and cds in hit_cds:
                    og_hits[og].add(gca)
                if cds.strip() != '':
                    og_total[og].add(gca)

for og in og_hits:
    print(og + '\t' + str(len(og_hits[og])) + '\t' + str(len(og_hits[og])/len(og_total[og])))
