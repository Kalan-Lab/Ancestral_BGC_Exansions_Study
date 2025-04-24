import os
import sys
from collections import defaultdict

antismash_file = 'AntiSMASH_Results.Strict.txt'
antismash_relaxed_file = 'AntiSMASH_Results.txt'
kofam_file = 'KOFamScan_Summary.txt'
cazy_file = 'CAZy_Annotations.txt'

gca_kofams = defaultdict(set)
with open(kofam_file) as okf:
    for line in okf:
        line = line.strip()
        ls = line.split('\t')
        gca_kofams[ls[0]].add(ls[3])

gca_total_caz = defaultdict(lambda: '0')
gca_distinct_caz = defaultdict(lambda: '0')
with open(cazy_file) as ocf:
    for line in ocf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('.')[0]
        total_cazy = ls[1]
        distinct_cazy = ls[2]
        gca_total_caz[gca] = total_cazy
        gca_distinct_caz[gca] = distinct_cazy

relaxed_bgcome_sizes = {}
relaxed_bgc_counts = {}
with open(antismash_relaxed_file) as oarf:
    for i, line in enumerate(oarf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        relaxed_bgcome_sizes[ls[0].split('.')[0]] = ls[4]
        relaxed_bgc_counts[ls[0].split('.')[0]] = ls[3]
with open(antismash_file) as oaf:
    for i, line in enumerate(oaf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0:
            print('\t'.join(ls + ['Relaxed_BGCome_Size', 'Relaxed_BGC_Count', 'Distinct_KOfams', 'Distinct_CAZy', 'Total_CAZy']))
        else:
            gca = ls[0].split('.')[0]
            print('\t'.join(ls + [str(relaxed_bgcome_sizes[gca]), str(relaxed_bgc_counts[gca]), str(len(gca_kofams[gca])), gca_distinct_caz[gca], gca_total_caz[gca]]))

