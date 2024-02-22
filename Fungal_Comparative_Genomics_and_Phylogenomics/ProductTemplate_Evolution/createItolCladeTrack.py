import os
import sys
from collections import defaultdict
from Bio import SeqIO

print('DATASET_COLORSTRIP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tClade')
print('DATA')
colors = {'Pezizomycotina': '#1f53a6', 'Other': '#000000'}

with open('PT_Seqs.faa') as opf:
    for rec in SeqIO.parse(opf, 'fasta'):
        clade = rec.id.split('|')[-1]
        if clade in colors:
            print(rec.id + '\t' + colors[clade] + '\t' + clade)
