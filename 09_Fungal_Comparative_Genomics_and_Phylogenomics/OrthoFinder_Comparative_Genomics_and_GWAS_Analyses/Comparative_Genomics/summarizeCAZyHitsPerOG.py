import os
import sys
from collections import defaultdict
from operator import itemgetter

evalue_threshold = 1e-17
coverage_threshold = 0.45

cazy_res_file = '../../Other_Annotations/CAZy/Results.txt'
og_file = 'Orthogroups.tsv'

data = []
with open(cazy_res_file) as ocrf:
    for line in ocrf:
        line = line.strip()
        ls = line.split('\t')
        cov = float(ls[6])
        evalue = float(ls[4])
        if evalue < evalue_threshold and cov > coverage_threshold:
            data.append([ls[0].split('.hmm')[0], ls[2], float(ls[5])])

top_hits = {}
for hit in sorted(data, key=itemgetter(2), reverse=True):
    if hit[1] in top_hits: continue
    top_hits[hit[1]] = hit[0]

with open(og_file) as oof:
    for i, line in enumerate(oof):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: continue
        og = ls[0]
        cazy = defaultdict(int)
        print_flag = False
        for lts in ls[1:]:
            for lt in lts.split(', '):
                if lt.strip() == '': continue
                if lt in top_hits:
                    cazy[top_hits[lt]] += 1
                    print_flag = True
                else:
                    cazy['no_cazy_match'] += 1
        
        if print_flag:
            cazy_string = []
            for j, fam in enumerate(sorted(cazy.items(), key=itemgetter(1), reverse=True)):
                cazy_string.append(fam[0] + '=' + str(fam[1]))
        
            print(og + '\t' + '; '.join(cazy_string))
