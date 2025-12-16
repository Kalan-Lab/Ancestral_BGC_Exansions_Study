import os
import sys
from Bio import SeqIO

all_prot_file = 'Proteins_with_HI_Assoc_Doms_GA.faa'
hmmsearch_file = 'HI_Full_Hits.Domain.Cut_GA.txt'

het_hits = set([])
with open(hmmsearch_file) as otcf:
    for line in otcf:
        if line.startswith('#'): continue
        line = line.strip()
        ls = line.split()
        if ls[3] == 'HET':
            het_hits.add(ls[0])

with open(all_prot_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        if rec.id in het_hits:
            print('>' + rec.description + '\n' + str(rec.seq))
