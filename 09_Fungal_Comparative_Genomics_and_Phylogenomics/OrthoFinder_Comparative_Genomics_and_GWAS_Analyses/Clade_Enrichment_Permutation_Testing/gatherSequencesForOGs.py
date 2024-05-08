import os
import sys
from Bio import SeqIO

sig_hogs = set([])
with open('Significant_Coarse_OGs.txt') as ohf:
    for i, line in enumerate(ohf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        sig_hogs.add(ls[1])

gcas = set([])
with open('Other_Pez.GCAs.txt') as of:
    for line in of:
        line = line.strip()
        gcas.add(line)

with open('BGC_Enriched_Clade.GCAs.txt') as of:
    for line in of:
        line = line.strip()
        gcas.add(line) 

cds_to_hog = {}
with open('../OrthoFinder_Results/Results_Apr09/Orthogroups/Orthogroups.tsv') as ohf:
    for i, line in enumerate(ohf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        hog = ls[0]
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

