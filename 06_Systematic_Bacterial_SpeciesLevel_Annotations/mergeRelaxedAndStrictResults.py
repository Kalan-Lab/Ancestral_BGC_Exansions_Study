import os
import sys
from collections import defaultdict

"""
0       Assembly_ID
1       Phylum
2       Genus
3       BGC_Count
4       Complete_BGC_Count
5       Genome_Size
6       BGC_Genome_Prop
7       Metallophore_Genome_Prop
8       NRPS_Genome_Prop
9       PKS_Genome_Prop
10      NRPS_or_PKS_Genome_Prop
11      NRPS_or_PKS_Prop
12      NRPS_or_PKS
13      Complete_NRPS_or_PKS
"""

relaxed_file = 'AntiSMASH_Stats.Updated_with_Complete_Counts.txt'
strict_file = 'AntiSMASH_Stats.Strict.txt'

gca_relaxed_data = {}
with open(relaxed_file) as orf:
    for i, line in enumerate(orf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca_relaxed_data[ls[0]] = ls[3:11] + ls[12:]

with open(strict_file) as osf:
    for i, line in enumerate(osf):
        line = line.strip()
        ls = line.split('\t') 
        if i == 0:
            print('\t'.join(['BGC_Sum_Relaxed', 'BGC_Sum_Strict', 'NRPS_or_PKS_Sum_Relaxed', 'NRPS_or_PKS_Sum_Strict'] + [x + '_Relaxed' for x in (ls[3:10] + ls[11:])] + [x + '_Strict' for x in (ls[3:10] + ls[11:])]))
        else:
            gca = ls[0]
            bss = str(float(ls[5])*float(ls[6]))
            bsr = str(float(gca_relaxed_data[gca][2])*float(gca_relaxed_data[gca][3]))
            npss = str(float(ls[5])*float(ls[10]))
            npsr = str(float(gca_relaxed_data[gca][2])*float(gca_relaxed_data[gca][7]))
            print('\t'.join([gca, bsr, bss, npsr, npss] + gca_relaxed_data[gca] + ls[3:11] + ls[12:]))
