import os
import sys
from collections import defaultdict
from operator import itemgetter

gtdb_r214_tsv = 'bac120_metadata_r214.tsv'

order_counts = defaultdict(int)
order_genera = defaultdict(set)
genus_counts = defaultdict(int)
with open(gtdb_r214_tsv) as ogrt:
    for i, line in enumerate(ogrt):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue 
        tax = ls[16]
        if not 'd__Bacteria;' in tax: continue
        order = tax.split(';o__')[1].split(';f__')[0]
        genus = tax.split(';g__')[1].split(';s__')[0]
        order_counts[order] += 1
        order_genera[order].add(genus)
        genus_counts[genus] += 1

top_order_genera = defaultdict(set)
for o in order_genera:
    top_genera = set([])
    top_genus_count = 0
    for g in order_genera[o]:
        gc = genus_counts[g]
        if gc > top_genus_count:
            top_genera = set([g])
            top_genus_count = gc
        elif gc == top_genus_count:
            top_genera.add(g)

    top_order_genera[o] = top_genera

assembly_level_values = {'Complete Genome': 4, 'Chromosome': 3, 'Scaffold': 2, 'Contig': 1}
date_informative_orders = set(['Vampirovibrionales', 'Cyanobacteriales'])

best_assemblies = defaultdict(lambda: [0, 0.0, None])
with open(gtdb_r214_tsv) as ogrt:
    for i, line in enumerate(ogrt):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        tax = ls[16]
        order = tax.split(';o__')[1].split(';f__')[0]
        if order_counts[order] >= 200 or order in date_informative_orders:
            gca = ls[54]
            n50 = int(ls[43])
            assembly_level = ls[45]
            assembly_level_val = assembly_level_values[assembly_level]
            if ls[43] == 'none': continue
            genus  = tax.split(';g__')[1].split(';s__')[0]
            if not genus in top_order_genera[order]: continue
            if not order in date_informative_orders:
                if n50 < 100000: continue
                if assembly_level_val != 4 and assembly_level_val != 3: continue
                if gca == 'GCA_024760545.1' or gca == 'GCA_024760525.1': continue
                if assembly_level_val > best_assemblies[order][0]:
                    best_assemblies[order] = [assembly_level_val, n50, gca]
                elif assembly_level_val == best_assemblies[order][0]:
                    if n50 > best_assemblies[order][1]:
                        best_assemblies[order] = [assembly_level_val, n50, gca]
            else:
                if assembly_level_val > best_assemblies[order][0]:
                    best_assemblies[order] = [assembly_level_val, n50, gca]
                elif assembly_level_val == best_assemblies[order][0]:
                    if n50 > best_assemblies[order][1]:
                        best_assemblies[order] = [assembly_level_val, n50, gca]

for order in sorted(order_counts.items(), key=itemgetter(1), reverse=True):
    if 'GCA' in '\t'.join([str(x) for x in best_assemblies[order[0]]]):
        print(order[0] + '\t' + str(order[1]) + '\t' + '\t'.join([str(x) for x in best_assemblies[order[0]]]))
