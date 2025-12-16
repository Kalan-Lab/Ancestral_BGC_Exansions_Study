import os
import sys
from Bio import SeqIO
from collections import defaultdict

all_prot_file = 'Proteins_with_HI_Assoc_Doms_GA.faa'
tc_cutoff_file = 'HI_Full_Hits.Domain.Cut_GA.txt'

coords = defaultdict(list)
with open(tc_cutoff_file) as otcf:
    for line in otcf:
        if line.startswith('#'): continue
        line = line.strip()
        ls = line.split()
        if ls[3] == 'HET':
            start = int(ls[17])
            end = int(ls[18])
            coords[ls[0]].append([start, end])

with open(all_prot_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        if rec.id in coords:
            for i, c in enumerate(coords[rec.id]):
                s, e = c
                print('>' + rec.description + '|' + str(s) + '|' + str(e) + '\n' + str(rec.seq)[s-1:e])
