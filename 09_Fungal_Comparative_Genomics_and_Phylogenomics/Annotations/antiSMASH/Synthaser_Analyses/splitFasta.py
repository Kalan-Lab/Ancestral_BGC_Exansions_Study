import os
import sys
from Bio import SeqIO

input_file = sys.argv[1]
prefix = sys.argv[2]

tmp = []
count = 1
with open(input_file) as of:
    for rec in SeqIO.parse(of, 'fasta'):
        if len(tmp) == 1000:
            outf = open(prefix + '.' + str(count) + '.faa', 'w')
            outf.write('\n'.join(tmp) + '\n')
            outf.close()
            tmp = []
            count += 1
        tmp.append('>' + rec.id)
        tmp.append(str(rec.seq))
if len(tmp) > 0:
    outf = open(prefix + '.' + str(count) + '.faa', 'w')
    outf.write('\n'.join(tmp) + '\n')
    outf.close()
