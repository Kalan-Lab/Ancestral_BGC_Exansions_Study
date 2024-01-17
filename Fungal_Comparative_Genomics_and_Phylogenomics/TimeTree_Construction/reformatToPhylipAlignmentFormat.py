import os
import sys
from Bio import SeqIO

print('156 13434')
for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    spaces = ' '*(42 - len(rec.id))
    print(rec.id + spaces + str(rec.seq))

