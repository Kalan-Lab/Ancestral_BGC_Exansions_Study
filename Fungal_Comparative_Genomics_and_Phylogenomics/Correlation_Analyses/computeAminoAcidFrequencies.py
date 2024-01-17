import os
import sys
from collections import defaultdict
from Bio import SeqIO

proteome_file = 'All_Fungal_Proteins.faa'

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

aa_counts = defaultdict(lambda: defaultdict(int))
with open(proteome_file) as opf:
    for rec in SeqIO.parse(opf, 'fasta'):
        gca = rec.id.split('|')[0]
        for r in str(rec.seq).upper():
            if r in aas:
                aa_counts[gca][r] += 1

print('sample\t' + '\t'.join(aas))
for gca in aa_counts:
    printlist = [gca]
    for aa in aa_counts[gca]:
        printlist.append(str(aa_counts[gca][aa]/float(sum(aa_counts[gca].values()))))
    print('\t'.join(printlist))
