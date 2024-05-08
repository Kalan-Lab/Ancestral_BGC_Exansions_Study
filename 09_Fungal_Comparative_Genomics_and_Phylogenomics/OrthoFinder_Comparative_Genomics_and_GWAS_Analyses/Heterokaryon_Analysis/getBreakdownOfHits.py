import os
import sys
from collections import defaultdict

gca_top_hits = defaultdict(lambda: [100000.0, set([])])
with open('Heterokaryon_Search_All_Prots.txt') as of:
    for line in of:
        line = line.strip()
        if line.startswith('#'): continue
        ls = line.split()
        evalue = float(ls[4])
        if evalue < 1e-3:
            cds = ls[0]
            hit = ls[3] 
            if gca_top_hits[cds][0] > evalue:
                gca_top_hits[cds] = [evalue, set([hit])]
            elif gca_top_hits[cds][0] == evalue:
                gca_top_hits[cds][1].add(hit)

for cds in gca_top_hits:
    print(cds + '\t' + str(gca_top_hits[cds][0]) + '\t' + str('; '.join(gca_top_hits[cds][1])))

