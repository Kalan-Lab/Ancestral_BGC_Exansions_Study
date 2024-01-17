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
        gc = 0
        tot = 0
        with open(genome) as ogf:
            for rec in SeqIO.parse(ogf, 'fasta'):
               for bp in str(rec.seq).upper():
                   if bp in set(['A', 'C', 'G', 'T']):
                       tot += 1
                       if bp in set(['C', 'G']):
                           gc += 1 
        gc_prop = gc/float(tot)
        print(ls[0] + '\t' + str(gc_prop))
