import os
import sys
from collections import defaultdict
from operator import itemgetter

gtdb_r214_tsv = 'bac120_metadata_r214.tsv' # can be downloaded from GTDB FTP
focal_phylum = sys.argv[1] # e.g. "p__Actinomycetota;"

genus_counts = defaultdict(int)
with open(gtdb_r214_tsv) as ogrt:
    for i, line in enumerate(ogrt):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue 
        tax = ls[16]
        if not focal_phylum in tax: continue
        genus = tax.split(';g__')[1].split(';s__')[0]
        genus_counts[genus] += 1

assembly_level_values = {'Complete Genome': 4, 'Chromosome': 3, 'Scaffold': 2, 'Contig': 1}

best_assemblies = defaultdict(lambda: [0, 0.0, None])
with open(gtdb_r214_tsv) as ogrt:
    for i, line in enumerate(ogrt):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        tax = ls[16]
        genus = tax.split(';g__')[1].split(';s__')[0]
        if genus_counts[genus] >= 5:
            assembly_level = ls[45]
            assembly_level_val = assembly_level_values[assembly_level]
            if ls[43] == 'none': continue
            n50 = int(ls[43])
            gca = ls[54]
            if assembly_level_val > best_assemblies[genus][0]:
                best_assemblies[genus] = [assembly_level_val, n50, gca]
            elif assembly_level_val == best_assemblies[genus][0]:
                if n50 > best_assemblies[genus][1]:
                    best_assemblies[genus] = [assembly_level_val, n50, gca]

for gen in sorted(genus_counts.items(), key=itemgetter(1), reverse=True):
    if gen[1] >= 5:
        print(gen[0] + '\t' + str(gen[1]) + '\t' + '\t'.join([str(x) for x in best_assemblies[gen[0]]]))
