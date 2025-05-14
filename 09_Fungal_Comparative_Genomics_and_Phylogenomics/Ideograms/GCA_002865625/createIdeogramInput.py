import os
import sys
from Bio import SeqIO
from operator import itemgetter
from collections import defaultdict

# Chrom    Start      End  Stain Arm

chrom_lengths = {}
with open('AntiSMASH_results/GCA_002865625/GCA_002865625.1_Rhier1_genomic.gbk') as of:
    for rec in SeqIO.parse(of, 'genbank'):
        chrom_lengths[rec.id] = len(str(rec.seq))

chrom_bgc_locs = defaultdict(list)
with open('BGC_Coordinates.txt') as obcf:
    for i, line in enumerate(obcf):
        line = line.strip()
        ls = line.split('\t')
        start = int(ls[2])
        end = int(ls[3])
        assert(start < end)
        chrom_bgc_locs[ls[1]].append([start, end])

print('\t'.join(['Chrom', 'Start', 'End', 'Stain', 'Arm']))
for chrom in chrom_lengths:
    cl = chrom_lengths[chrom]
    if len(chrom_bgc_locs[chrom]) == 0:
        print('\t'.join([str(x) for x in [chrom, 0, cl, 'background', '1']]))
        continue
    first_start = min([x[0] for x in chrom_bgc_locs[chrom]])
    if first_start > 0:
        print('\t'.join([str(x) for x in [chrom, 0, first_start, 'background', '1']]))
    
    last_end = None
    for i, bgc in enumerate(sorted(chrom_bgc_locs[chrom], key=itemgetter(1))):
        start = bgc[0]
        end = bgc[1]
        print('\t'.join([str(x) for x in [chrom, start, end, 'bgc', '1']]))
        if i+1 < len(chrom_bgc_locs[chrom]):
            next_start = sorted(chrom_bgc_locs[chrom], key=itemgetter(1))[i+1][0]
            print('\t'.join([str(x) for x in [chrom, end, next_start, 'background', '1']]))
        last_end = end
    if last_end != None and last_end <= cl:
        print('\t'.join([str(x) for x in [chrom, last_end, cl, 'background', '1']]))
