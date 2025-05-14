import os
import sys
from Bio import SeqIO

prot_locs = {}
with open('GCA_000827335.faa') as opf:
    for rec in SeqIO.parse(opf, 'fasta'):
        desc = rec.description
        prot = rec.id
        chrom = desc.split()[1]
        midpoint = ((int(desc.split()[2])-1) + (int(desc.split()[3])-1))/2.0
        prot_locs[prot] = [chrom, midpoint]

print('\t'.join(['Protein', 'Chrom', 'Position', 'Type']))
with open('HIT_Hits.Domain.Cut_TC.txt') as ohf:
    for line in ohf:
        if line.startswith('#'): continue
        line = line.strip()
        ls = line.split()
        hit_type = ls[3]
        prot = ls[0]
        chrom, pos = prot_locs[prot]
        print('\t'.join([prot, chrom, str(pos), hit_type]))
