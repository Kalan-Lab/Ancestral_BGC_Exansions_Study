import os
import sys
from Bio import SeqIO

all_prot_file = '../All_Proteins.faa'
tc_cutoff_file = '../HIT_Hits.Cut_TC.txt'

het_hits = set([])
with open(tc_cutoff_file) as otcf:
    for line in otcf:
        if line.startswith('#'): continue
        line = line.strip()
        ls = line.split()
        if ls[2] == 'NACHT':
            het_hits.add(ls[0])

with open(all_prot_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        if rec.id in het_hits:
            print('>' + rec.description + '\n' + str(rec.seq))
