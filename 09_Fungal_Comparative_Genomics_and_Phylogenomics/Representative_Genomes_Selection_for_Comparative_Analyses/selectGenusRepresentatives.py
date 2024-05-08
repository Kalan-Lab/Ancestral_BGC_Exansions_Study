import os
import sys
from collections import defaultdict

as_file = 'AntiSMASH_Stats.Updated_with_Complete_Counts.txt'
gg_file = 'GCA_to_Genus.txt'

gca_to_bgcomesum = {}
with open(as_file) as oaf:
    for i, line in enumerate(oaf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca_to_bgcomesum[ls[0]] = float(ls[4])

genus_biggest_bgcome_sum = defaultdict(lambda: ['NA', -1.0])
with open(gg_file) as ogf:
    for line in ogf:
        line = line.strip()
        gca, genus = line.split('\t')
        bgcome_sum = gca_to_bgcomesum[gca]
        if genus_biggest_bgcome_sum[genus][1] < bgcome_sum:
            genus_biggest_bgcome_sum[genus] = [gca, bgcome_sum]

for genus in genus_biggest_bgcome_sum:
    print(genus_biggest_bgcome_sum[genus][0] + '\t' + genus)
