import os
import sys
from collections import defaultdict
from Bio import SeqIO

prot_dir = os.path.abspath('Proteomes/') + '/'
antismash_stats_file = 'AntiSMASH_Stats.txt'
cazy_stats_file = 'CAZy_Summary.txt'
het_stats_file = 'HI.iTol.txt'
clade_files = ['Clades/Neocallimastigomycota.txt', 'Clades/Agaricomycetes.txt', 'Clades/BGC_Enriched_Pezizomycotina.txt']

name_to_clade = defaultdict(lambda: 'Other')
for cf in clade_files:
    clade = cf.split('.')[0].split('/')[-1]
    with open(cf) as ocf:
        for line in ocf:
            line = line.strip()
            name_to_clade[line] = clade

proteins = defaultdict(set)
for f in os.listdir(prot_dir):
    with open(prot_dir + f) as opf:
        for rec in SeqIO.parse(opf, 'fasta'):
            proteins[f.split('.faa')[0]].add(rec.id)

cazy = defaultdict(set)
with open(cazy_stats_file) as ocsf:
    for line in ocsf:
        line = line.strip()
        ls = line.split('\t')
        cazy[ls[0]] = ls[1]

hi_counts = defaultdict(int)
with open(het_stats_file) as ohsf:
    for i, line in enumerate(ohsf):
        if i <= 5: continue
        line = line.strip()
        ls = line.split('\t')
        hi_counts[ls[0]] = int(ls[1]) + int(ls[2])

print('\t'.join(['name', 'gca', 'clade', 'bgcome_size', 'genome_size', 'cazy_count', 'hi_count', 'protein_count']))
with open(antismash_stats_file) as oasf:
    for line in oasf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        name = ls[1]
        bs = ls[4]
        gs = ls[5]
        print('\t'.join([name, gca, name_to_clade[name], bs, gs, str(cazy[gca]), str(hi_counts[name]), str(len(proteins[gca]))]))

