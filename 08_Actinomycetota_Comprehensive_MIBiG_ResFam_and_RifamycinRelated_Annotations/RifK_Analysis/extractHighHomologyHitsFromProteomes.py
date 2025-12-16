import os
import sys
from Bio import SeqIO

search_dir = 'RifK_Search/'
proteomes_dir = 'Proteomes/'

for f in os.listdir(search_dir):
    ko_file = search_dir + f
    prot_file = proteomes_dir + f.split('.txt')[0]

    matches = {}
    hit = set([])
    with open(ko_file) as okf:
        for i, line in enumerate(okf):
            if line.startswith('#'): continue
            line = line.strip()
            ls = line.split()
            if line.startswith('*'):
                matches[ls[1]] = ls[2]
                hit.add(ls[1])
            elif float(ls[4]) < 1e-100: # original was 1e-10
                matches[ls[0]] = ls[1]

    gca = f.split('.')[0]
    with open(prot_file) as opf:
        for rec in SeqIO.parse(opf, 'fasta'):
            if rec.id in matches:
                flag = ''
                if rec.id in hit:
                    flag = 'hit'
                print('>' + gca + '|' + rec.id + '|' + matches[rec.id] + '|' + flag + '\n' + str(rec.seq))
