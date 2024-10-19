import os
import sys
from collections import defaultdict
from Bio import SeqIO

select = set([])
with open('Terpene_OGs_in_atl25perGenomes.txt') as of:
    for line in of:
        line = line.strip()
        ls = line.split('\t')
        select.add(ls[0])

cds = set([])
cds_to_og = {}
with open('Orthogroups.tsv') as ogf:
    for i, line in enumerate(ogf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        og = ls[0]
        if og in select:
            for cs in ls[1:]:
                for c in cs.split(', '):
                    if c.strip() != '':
                        cds.add(c)
                        cds_to_og[c] = og

prot_dir = 'Proteomes/'
for f in os.listdir(prot_dir):
    with open(prot_dir + f) as opf:
        for rec in SeqIO.parse(opf, 'fasta'):
            if rec.id in cds:
                print('>' + cds_to_og[rec.id] + '|' + rec.id + '\n' + str(rec.seq))
