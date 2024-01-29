import os
import sys
from collections import defaultdict

aafreqs_file = 'Fungal_Amino_Acid_Frequencies.txt'
gc_file = 'Fungal_Genome_GC_Contents.txt'
gs_file = 'Fungal_Genomes.txt'
enriched_clade_file = 'BGC_Enriched_Clade_Genomes.txt'
ascomycota_file = 'Ascomycota_Genomes.txt'
as_stats_file = 'AntiSMASH_Stats.Updated.txt'
starship_file = 'Starship.iTol.txt'
cazy_file = 'CAZy.Distinct.iTol.txt'

cazy_counts = defaultdict(int)
with open(cazy_file) as ocf:
    for i, line in enumerate(ocf):
        line = line.strip()
        if i > 5:
            ls = line.split('\t')
            gca = '_'.join(ls[0].split('_')[-2:])
            cazy_counts[gca] = int(ls[1])

starship_counts = defaultdict(int)
with open(starship_file) as ocf:
    for i, line in enumerate(ocf):
        line = line.strip()
        if i > 5:
            ls = line.split('\t')
            gca = '_'.join(ls[0].split('_')[-2:])
            starship_counts[gca] = int(ls[1])

ascomycota = set([])
enriched_clade = set([])

with open(ascomycota_file) as oaf:
    for line in oaf:
        line = line.strip()
        ascomycota.add('_'.join(line.split('_')[-2:]))

with open(enriched_clade_file) as oef:
    for line in oef:
        line = line.strip()
        enriched_clade.add('_'.join(line.split('_')[-2:]))

bgc_sums = defaultdict(int)
bgc_props = defaultdict(int)
genome_sizes = defaultdict(int)
with open(as_stats_file) as oasf:
    for i, line in enumerate(oasf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        bgc_prop = ls[5] 
        bgc_sum = ls[3]
        gca = ls[0]
        gs = ls[4]
        bgc_sums[gca] = bgc_sum
        bgc_props[gca] = bgc_prop
        genome_sizes[gca] = gs

aa_freqs = {}
with open(aafreqs_file) as oaf:
    for i, line in enumerate(oaf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        aa_freqs[ls[0]] = ls[1:]

gcs = {}
with open(gc_file) as ogf:
    for line in ogf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        gcs[gca] = ls[1]

print('\t'.join(['gca', 'clade', 'bgc_sum', 'bgc_prop', 'genome_size', 'gc', 'cazy_counts', 'starship_counts', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])) 

with open(gs_file) as ogf:
    for line in ogf:
        gca = line.strip()
        genome_size = genome_sizes[gca]
        is_ascomycota = 'Non-Ascomycota'
        if gca in ascomycota:
            is_ascomycota = 'Ascomycota'
            if gca in enriched_clade:
                is_ascomycota = 'Ascomycota - BGC Enriched Clade'
        bgc_sum = bgc_sums[gca]
        bgc_prop = bgc_props[gca]
        cazy_count = cazy_counts[gca]
        starship_count = starship_counts[gca]
        gc = gcs[gca]
        #if is_ascomycota == 'Ascomycota - BGC Enriched Clade':
        print('\t'.join([str(x) for x in [gca, is_ascomycota, bgc_sum, bgc_prop, genome_size, gc, cazy_count, starship_count] + aa_freqs[gca]]))
