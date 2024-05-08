import os
import sys
from collections import defaultdict
from Bio import SeqIO

rae_hit_dir = 'BLASTn_Results/'
gca_to_tax_file = '../GCA_to_Taxonomy.txt'
tree_fam_file = 'FamilyTree_Represented.txt'
rae_fna_file = 'RAE_Palindrome.fna'

gca_to_fam = {}
fam_to_cla = {}
with open(gca_to_tax_file) as ogttf:
    for line in ogttf:
        line = line.strip()
        ls = line.split('\t')
        fam = ls[1].split('f__')[1].split(';g__')[0]
        cla = ls[1].split('c__')[1].split(';o__')[0]
        gca_to_fam[ls[0]] = fam
        fam_to_cla[fam] = cla

id_to_seq = {}
with open(rae_fna_file) as orff:
    for rec in SeqIO.parse(orff, 'fasta'):
        id_to_seq[rec.id] = str(rec.seq)

fam_hits = defaultdict(int)
fam_total = defaultdict(int)
for f in os.listdir(rae_hit_dir):
    rae_hit_file = rae_hit_dir + f
    hit = 0
    seqs_found = set([])
    with open(rae_hit_file) as orhf:
        for line in orhf:
            line = line.strip()
            ls = line.split('\t')
            if float(ls[2]) == 100.0 and float(ls[-1]) == 100.0:
                hit += 1
                seqs_found.add(id_to_seq[ls[0]])
    gca = '_'.join(f.split('_')[:2]).split('.')[0]
    if len(seqs_found) > 0:
        print(gca + '\t' + ', '.join(seqs_found))
