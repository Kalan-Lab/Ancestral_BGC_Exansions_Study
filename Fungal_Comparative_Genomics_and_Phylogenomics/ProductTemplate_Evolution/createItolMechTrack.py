import os
import sys
from collections import defaultdict
from Bio import SeqIO

print('DATASET_PIECHART')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tClade')
print('FIELD_LABELS\tC2-C7\tC4-C9\tOther')
print('FIELD_COLORS\t#edd6a8\t#614309\t#787878')
print('DATA')

with open('PT_Seqs.faa') as opf:
    for rec in SeqIO.parse(opf, 'fasta'):
        if rec.id.startswith('Ref'):
            mech = rec.id.split('|')[2]
            printlist = [rec.id, '0', '1']
            if mech == 'C2-C7':
                printlist += ['1', '0', '0']
            elif mech == 'C4-C9':
                printlist += ['0', '1', '0']
            else:
                printlist += ['0', '0', '1']

            print('\t'.join(printlist))
