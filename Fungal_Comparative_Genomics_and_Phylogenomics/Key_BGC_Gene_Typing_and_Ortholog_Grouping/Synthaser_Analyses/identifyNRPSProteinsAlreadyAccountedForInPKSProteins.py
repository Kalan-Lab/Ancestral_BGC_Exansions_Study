import os
import sys
from Bio import SeqIO

pks_prots = set([])
with open('../All_PKS_Full_Proteins.faa') as opf:
    for rec in SeqIO.parse(opf, 'fasta'):
        pks_prots.add(' '.join(rec.description.split()[1:]))

with open('../All_NRPS_Full_Proteins.faa') as onf:
    for rec in SeqIO.parse(onf, 'fasta'):
        prot_id = ' '.join(rec.description.split()[1:])
        if prot_id in pks_prots:
            print(rec.description)
