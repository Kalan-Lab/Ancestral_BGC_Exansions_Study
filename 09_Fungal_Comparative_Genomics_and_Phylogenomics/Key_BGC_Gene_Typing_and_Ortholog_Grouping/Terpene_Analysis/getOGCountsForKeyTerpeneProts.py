import os
import sys
from Bio import SeqIO
from collections import defaultdict

outgroups = set(['GCA_000151315', 'GCA_000002865', 'GCA_000001215', 'GCA_000388065'])

gca_seqs = defaultdict(set)
with open('All_Key_Terpene_Proteins.faa') as oktf:
    for rec in SeqIO.parse(oktf, 'fasta'):
        rec_desc = rec.description
        gca = rec_desc.split()[1].split('|')[1]
        gca_seqs[gca].add(str(rec.seq))

prot_dir = 'Proteomes/'
cdss = set([])
for f in os.listdir(prot_dir):
    gca = f.split('.faa')[0]
    with open(prot_dir + f) as opf:
        for rec in SeqIO.parse(opf, 'fasta'):
            if rec.seq in gca_seqs[gca]:
                cdss.add(rec.id)

hog_total_counts = defaultdict(set)
hog_terp_counts = defaultdict(set)
gcas = []
with open('Orthogroups.tsv') as ohf:
    for i, line in enumerate(ohf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: 
            gcas = ls[1:]
            continue
        hog = ls[0]#.split('N0.')[1]
        for j, cds in enumerate(ls[1:]):
            gca = gcas[j]
            if gca in outgroups: continue
            for c in cds.split(', '):
                if c.strip() != '':
                    hog_total_counts[hog].add(gca)
                    if c in cdss:
                        hog_terp_counts[hog].add(gca)

for hog in hog_terp_counts:
    if len(hog_terp_counts[hog]) >= 63:
        print(hog + '\t' + str(len(hog_terp_counts[hog])/len(hog_total_counts[hog])) + '\t' + str(len(hog_terp_counts[hog])) + '\t' + str(len(hog_total_counts[hog])))
