import os
import sys
from Bio import SeqIO

overview_file = 'Overview_File.txt'

with open(overview_file) as ovf:
    for i, line in enumerate(ovf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        genome = ls[-2]
        tot = 0
        with open(genome) as ogf:
            for rec in SeqIO.parse(ogf, 'fasta'):
                tot += len(rec.seq)
        print(ls[0] + '\t' + str(tot))
