import os
import sys
from collections import defaultdict
from Bio import SeqIO

gca_to_class = {}
gca_to_fam = {}
with open('GCA_to_Taxonomy.txt') as ogtf:
    for line in ogtf:
        line = line.strip()
        gca, tax = line.split('\t')
        tax_class = tax.split(';c__')[1].split(';')[0]
        tax_fam = tax.split(';f__')[1].split(';')[0]
        if tax_class == 'Actinomycetia':
            gca_to_class[gca] = tax_class
        else:
            gca_to_class[gca] = 'Other Actinomycetota class'
        gca_to_fam[gca] = tax_fam


top_hits = defaultdict(lambda: [set([]), 100000.0])
with open("DAHP_Typing_Full.txt") as ott:
    for line in ott:
        line = line.strip()
        ls = line.split()
        if line.startswith('#'): continue
        evalue = float(ls[4])
        hit = ls[2]
        prot = ls[0]
        if top_hits[prot][1] > evalue:
            top_hits[prot] = [set([hit]), evalue]
        elif top_hits[prot][1] == evalue:
            top_hits[prot][0].add(hit)

print('Enzyme\tClass\tCount')

dhap_faa_file = 'DAHP_Synthase.faa'
rifh_fai_file = 'Filtered_Fai_Hits_with_RifH.txt'

clade_class_counts = defaultdict(lambda: defaultdict(set))

with open(dhap_faa_file) as odff:
    for rec in SeqIO.parse(odff, 'fasta'):
        gca = rec.id.split('|')[0].split('.')[0]
        gca_class = gca_to_class[gca]        
        gca_class_group = 'Actinomycetia'
        if not gca_class == 'Actinomycetia':
            gca_class_group = 'Other Actinomycetota class'
        dtype = 'Not matching either type 1 or 2 DAHP synthase TIGR profiles'
        if rec.id in top_hits:
            hits = top_hits[rec.id][0]
            assert(len(hits) == 1)
            dtype = list(hits)[0]
        clade_class_counts[dtype][gca_class_group].add(gca)

with open(rifh_fai_file) as orff:
    for line in orff:
        line = line.strip()
        gca = line.split('\t')[0].split('.')[0]
        gca_class = gca_to_class[gca]
        gca_class_group = 'Actinomycetia'
        if not gca_class == 'Actinomycetia':
            gca_class_group = 'Other Actinomycetota class'
        clade_class_counts['RifH homolog from AHBA neighborhood'][gca_class_group].add(gca)

print(clade_class_counts.keys())
for clade in clade_class_counts:
    for clas in ['Actinomycetia', 'Other Actinomycetota class']:
        tot_class_count = 0
        if clas == 'Actinomycetia':
            tot_class_count = 28281
        else:
            tot_class_count = 5523
        print(clade + '\t' + clas + '\t' + str(len(clade_class_counts[clade][clas])/float(tot_class_count)))
