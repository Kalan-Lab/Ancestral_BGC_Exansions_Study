import os
import sys
from collections import defaultdict
from Bio import SeqIO

protein_to_info = defaultdict(lambda: defaultdict(list))
with open('All_BGC_Proteins.txt') as of:
    for line in of:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        prot = ls[-1]
        prot_info = ls[:-1]
        protein_to_info[gca][prot].append(prot_info)

cds_to_og = {}
with open('OrthoFinder_Results/Results_Apr09/Phylogenetic_Hierarchical_Orthogroups/N0.tsv') as onof:
    for i, line in enumerate(onof):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        for gs in ls[3:]:
            for g in gs.split(', '):
                if g != '':
                    cds_to_og[g] = ls[0]

prot_dir = 'Proteomes/'
for f in os.listdir(prot_dir):
    gca = f.split('.')[0]
    with open(prot_dir + f) as opf:
        for rec in SeqIO.parse(opf, 'fasta'):
            rid = rec.id
            rseq = str(rec.seq)
            if rseq in protein_to_info[gca]:
                val = protein_to_info[gca][rseq]
                og = 'NA'
                if rid in cds_to_og:
                    og = cds_to_og[rid]
                is_uniq = True
                if len(val) > 1: is_uniq = False
                print(og + '\t' + rid +'\t' + str(is_uniq) + '\t' + '\t'.join(['\t'.join(x) for x in val]))
