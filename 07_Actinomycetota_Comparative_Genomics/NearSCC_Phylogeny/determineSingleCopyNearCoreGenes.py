import os
import sys

with open('Orthogroups.tsv') as otf:
    for i, line in enumerate(otf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: continue
        og = ls[0]
        vals = []
        for j, gs in enumerate(ls[1:]):
            gene_count = gs.split(', ')
            if gs == '': vals.append(0)
            else: vals.append(len(gene_count))
        if max(vals) > 1: continue
        prop_core = sum(vals)/len(vals)
        if prop_core >= 0.8:
            print(og + '\t' + str(prop_core))
        
