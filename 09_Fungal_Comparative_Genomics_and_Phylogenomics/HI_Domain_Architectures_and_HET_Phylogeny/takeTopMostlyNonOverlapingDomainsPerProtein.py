import os
import sys
from collections import defaultdict
from operator import itemgetter

domain_annot_file = sys.argv[1]

protein_dom_hits = defaultdict(list)
with open(domain_annot_file) as odaf:
    for line in odaf:
        line = line.strip()
        prot, dom, start, end, score, _  = line.split('\t')
        protein_dom_hits[prot].append([dom, int(start), int(end), float(score)])

protein_dom_hits_retained = defaultdict(list)
for prot in protein_dom_hits:
    coords_hit_already = set([])
    for di in sorted(protein_dom_hits[prot], key=itemgetter(3), reverse=True):
        overlapping_better_hits = 0
        dom_coords = set([])
        for pp in range(di[1], di[2]+1):
            dom_coords.add(pp)
            if pp in coords_hit_already:
                overlapping_better_hits += 1
        prop_overlapping = overlapping_better_hits/len(dom_coords)
        if prop_overlapping < 0.1:
            protein_dom_hits_retained[prot].append(di)
            coords_hit_already = coords_hit_already.union(dom_coords)

for prot in protein_dom_hits_retained:
    domain_architecture = []
    for di in sorted(protein_dom_hits_retained[prot], key=itemgetter(1), reverse=False):
        if len(domain_architecture) == 0 or domain_architecture[-1] != di[0]:
            domain_architecture.append(di[0])
    print(prot + '\t' + '|'.join(domain_architecture))
        
