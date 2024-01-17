import os
import sys
from Bio import SeqIO
from collections import defaultdict

faa_msa_dir = 'Core_Protein_MSAs_Trimmed/'

gca_to_species_name = {}
with open('Overview_File.txt') as ovf:
    for line in ovf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        species_name_simp = '_'.join(ls[1].split(';')[-1].split()).split('_')[0]
        gca_to_species_name[gca] = (species_name_simp  + '_'+ gca).replace('(', '_').replace(')', '_')


sample_msa_seqs = defaultdict(lambda: '')
for f in os.listdir(faa_msa_dir):
    with open(faa_msa_dir + f) as off:
        for rec in SeqIO.parse(off, 'fasta'):
            sample_msa_seqs[rec.id] += str(rec.seq) 

for g in sample_msa_seqs:
    print('>' + gca_to_species_name[g] + '\n' + str(sample_msa_seqs[g]))
