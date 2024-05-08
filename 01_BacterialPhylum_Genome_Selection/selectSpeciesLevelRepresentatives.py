import os
import sys
from collections import defaultdict
from operator import itemgetter

gtdb_r214_tsv = 'bac120_metadata_r214.tsv' # can be downloaded from GTDB FTP
genus_rep_file = sys.argv[1] # listing of genus representatives for phylum from which to select species-level representatives
focal_phylum = sys.argv[2] # identifier for focal phylum, e.g. p__Actinomycetota; 

genera = set([])
with open(genus_rep_file) as oif:
    for line in oif:
        line = line.strip()
        ls = line.split('\t')
        genera.add(ls[0])

bad_gcas = set([])
try:
    """
    File with problematic GCAs for assemblies that were unable to be downloaded from NCBI -
    likely because assembly authors made it private retroactively following incorporation 
    of the assembly into GTDB.
    """
    for line in open('Bad_GCAs.txt'):
        line = line.strip()
        bad_gcas.add(line)
except:
    pass

assembly_level_values = {'Complete Genome': 4, 'Chromosome': 3, 'Scaffold': 2, 'Contig': 1}
genus_species_assemblies = defaultdict(lambda: defaultdict(list))
with open(gtdb_r214_tsv) as ogrt:
    for i, line in enumerate(ogrt):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue 
        tax = ls[16]
        genus = tax.split(';g__')[1].split(';s__')[0]
        species = tax.split(';s__')[-1]
        if not genus in genera: continue
        if ls[43] == 'none': continue
        n50 = int(ls[43])
        gca = ls[54]
        if n50 <= 100000: continue
        assembly_level = ls[45]
        if gca in bad_gcas: continue
        assembly_level_val = assembly_level_values[assembly_level]
        genus_species_assemblies[genus][species].append([assembly_level_val, n50, gca])

for g in genus_species_assemblies:
    for s in genus_species_assemblies[g]:
        for i, ga in enumerate(sorted(genus_species_assemblies[g][s], key = lambda x: (x[0], x[1]), reverse=True)):
            if i == 0:
                print(focal_phylum + '\t' + g + '\t' + s + '\t' + '\t'.join([str(x) for x in ga])) 
