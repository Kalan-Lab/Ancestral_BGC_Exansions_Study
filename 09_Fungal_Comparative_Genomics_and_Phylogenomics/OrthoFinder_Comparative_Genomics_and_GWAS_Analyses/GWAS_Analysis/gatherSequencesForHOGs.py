import os
import sys
from Bio import SeqIO

sig_hogs = set([])
with open('Significant_Associations_Dikarya.txt') as ohf:
    for i, line in enumerate(ohf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        sig_hogs.add(ls[0].split('N0.')[1])

gcas = set([])
with open('Dikarya_Names.txt') as of:
    for line in of:
        line = line.strip()
        gcas.add('_'.join(line.split('_')[-2:]))

cds_to_hog = {}
with open('Hierarchical_Orthogroups.tsv') as ohf:
    for i, line in enumerate(ohf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        hog = ls[0].split('N0.')[1]
        if not hog in sig_hogs: continue
        for gs in ls[1:]:
            for g in gs.split(', '):
                if g.strip() != '':
                    cds_to_hog[g] = hog

prot_dir = '../Proteomes/'
for f in os.listdir(prot_dir):
    gca = f.split('.faa')[0]
    if gca in gcas:
        for rec in SeqIO.parse(prot_dir + f, 'fasta'):
            if rec.id in cds_to_hog:
                print('>' + cds_to_hog[rec.id] + '|' + rec.id + '\n' + str(rec.seq))

