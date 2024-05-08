import os
import sys
from collections import defaultdict
from operator import itemgetter

gtdb_r214_tsv = sys.argv[1] # either 'bac120_metadata_r214.tsv' or 'ar53_metadata_r214.tsv' from GTDB.
individual_orders_file = sys.argv[2] # file listing selected orders (e.g. 74 in the final phylogeny).

orders = set([])
with open(individual_orders_file) as oif:
    for line in oif:
        line = line.strip()
        ls = line.split('\t')
        orders.add(ls[0])

bad_gcas = set([])
for line in open('Bad_GCAs.txt'):
    line = line.strip()
    bad_gcas.add(line)

assembly_level_values = {'Complete Genome': 4, 'Chromosome': 3, 'Scaffold': 2, 'Contig': 1}
order_species_assemblies = defaultdict(lambda: defaultdict(list))
with open(gtdb_r214_tsv) as ogrt:
    for i, line in enumerate(ogrt):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue 
        tax = ls[16]
        order = tax.split(';o__')[1].split(';f__')[0]
        if not order in orders: continue
        if ls[43] == 'none': continue
        n50 = int(ls[43])
        gca = ls[54]
        genus = tax.split(';g__')[1].split(';s__')[0]
        if n50 <= 100000: continue
        assembly_level = ls[45]
        if gca in bad_gcas: continue
        assembly_level_val = assembly_level_values[assembly_level]
        order_species_assemblies[order][genus].append([assembly_level_val, n50, gca])

for o in order_species_assemblies:
    for g in order_species_assemblies[o]:
        for i, ga in enumerate(sorted(order_species_assemblies[o][g], key = lambda x: (x[0], x[1]), reverse=True)):
            if i == 0:
                print(o + '\t' + g + '\t' + '\t'.join([str(x) for x in ga])) 
