import os
import sys
from Bio import SeqIO

seqs = []
for rec in SeqIO.parse('Core_Alignment_Concatenated.faa', 'fasta'):
    if rec.id == 'Ilumatobacter_A' or rec.id == 'AWTP1-35':
        seqs.append(list(str(rec.seq)))

for pos, als in enumerate(zip(*seqs)):
    als = list(als)
    if len(set(als)) == 1:
        if als[0] != '-':
            print(str(pos) + '\t' + str(als[0]))
