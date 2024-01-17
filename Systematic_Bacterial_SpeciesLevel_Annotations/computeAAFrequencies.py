import os
import sys
from Bio import SeqIO
from collections import defaultdict

aa_set = set(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'O', 'S' , 'U', 'T', 'W', 'Y', 'V'])

genome_aa_counts = defaultdict(lambda: defaultdict(int))
with open(sys.argv[1]) as otgf:
    for rec in SeqIO.parse(otgf, 'fasta'):
        genome = rec.id.split('.')[0]
        seq = str(rec.seq).upper()
        for aa in seq:
            if aa in aa_set:
                genome_aa_counts[genome][aa] += 1

print('genome\t' + '\t'.join(sorted(list(aa_set))))
for g in genome_aa_counts:
    gtot = sum(genome_aa_counts[g].values())
    row = [g]
    for aa in sorted(list(aa_set)):
        aac = genome_aa_counts[g][aa]
        aaf = aac/float(gtot)
        row.append(str(aaf))
    print('\t'.join(row))
