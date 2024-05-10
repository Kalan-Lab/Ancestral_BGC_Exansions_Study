import os
import sys
import numpy as np
from Bio import SeqIO

gaps = []
with open(sys.argv[1]) as of:
    for rec in SeqIO.parse(of, 'fasta'):
        gaps.append(len([1 for x in str(rec.seq) if x == '-' or x.upper() == 'X']))

#print('\n'.join([str(x) for x in gaps]))

#lower_25_count = np.percentile(sorted(gaps), 75)

gap_limit = 125

with open(sys.argv[1]) as of:
    for rec in SeqIO.parse(of, 'fasta'):
        gap_count = len([1 for x in str(rec.seq) if x == '-' or x.upper() == 'X'])
        if gap_count < gap_limit:
            print('>' + rec.description + '\n' + str(rec.seq))
