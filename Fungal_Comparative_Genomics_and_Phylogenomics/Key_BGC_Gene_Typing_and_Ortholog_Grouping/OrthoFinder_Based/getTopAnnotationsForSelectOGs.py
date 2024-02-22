import os
import sys
from collections import defaultdict

annot_file = 'DIAMOND_to_MIBiG.txt'
og_file = 'OrthoFinder_Results/Results_Feb12/Orthogroups/Orthogroups.tsv'
plot_file = 'Plotting_Input.txt'


select_ogs = set([])
with open(plot_file) as opf:
    for line in opf:
        line = line.strip()
        ls = line.split('\t')
        select_ogs.add(ls[0])

top_blast_hits = defaultdict(lambda: [0.0, '0.0', 'NA'])
with open(annot_file) as oaf:
    for line in oaf:
        line = line.strip()
        ls = line.split('\t')
        bs = float(ls[2])
        query = ls[0]
        hit = ls[1]
        scov = ls[-1]
        if bs > top_blast_hits[query][0]:
            top_blast_hits[query] = [bs, scov, hit]

with open(og_file) as oof:
    for i, line in enumerate(oof):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        og = ls[0]
        if not og in select_ogs: continue
        for ps in ls[1:]:
            for p in ps.split(', '):
                print(og + '\t' + p + '\t' +  str(top_blast_hits[p][0]) + '\t'+  top_blast_hits[p][1] + '\t' + top_blast_hits[p][2])
