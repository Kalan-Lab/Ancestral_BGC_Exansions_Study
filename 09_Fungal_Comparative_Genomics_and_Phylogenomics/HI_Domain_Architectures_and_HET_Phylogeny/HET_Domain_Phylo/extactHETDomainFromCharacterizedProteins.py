import os
import sys
from Bio import SeqIO
from collections import defaultdict

all_prot_file = 'Characterized_Proteins/Characterized_HI-Associated_Proteins.faa'
annot_file = 'Characterized_Proteins/Characterized_HI-Associated_Proteins.Pfam_Annotations_Gathering.txt'

coords = defaultdict(list)
with open(annot_file) as oaf:
    for line in oaf:
        if line.startswith('#'): continue
        line = line.strip()
        ls = line.split()
        if ls[1] == 'HET':
            start = int(ls[2])
            end = int(ls[3])
            coords[ls[0]].append([start, end])

with open(all_prot_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        if rec.id in coords:
            for i, c in enumerate(coords[rec.id]):
                s, e = c
                print('>' + rec.description + '|' + str(s) + '|' + str(e) + '\n' + str(rec.seq)[s-1:e])
