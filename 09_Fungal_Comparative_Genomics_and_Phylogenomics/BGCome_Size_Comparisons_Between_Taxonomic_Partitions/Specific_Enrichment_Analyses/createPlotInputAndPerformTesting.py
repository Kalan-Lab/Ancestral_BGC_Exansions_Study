import os
import sys
from scipy import stats
from collections import defaultdict

antismash_stats_file = 'AntiSMASH_Stats.Updated_WithoutContaminants.txt'

clade_1_file = sys.argv[1]
clade_2_file = sys.argv[2]
grid = sys.argv[3]

genus_reps_file = 'Genus_Representative_GCAs.txt'
genus_reps = set([])
with open(genus_reps_file) as ogrf:
    for line in ogrf:
        line = line.strip()
        genus_reps.add(line)


clade_1 = set([])
clade_2 = set([])
gca_to_clade = {}

with open(clade_1_file) as oc1f:
    for line in oc1f:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        if gca in genus_reps:
            clade_1.add(gca)
            gca_to_clade[gca] = 'clade_1'

with open(clade_2_file) as oc2f:
    for line in oc2f:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        if gca in genus_reps:
            clade_2.add(gca)
            gca_to_clade[gca] = 'clade_2'

clade_1_bgcomes = []
clade_2_bgcomes = []

print('GCA\tPhylum_or_Clade\tComplete_BGC_Count\tBGCome_Size\tStrict_NRPS_or_PKSome_Size\tGrid')
with open(antismash_stats_file) as oasf:
    for line in oasf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        if not gca in gca_to_clade: continue
        clade = gca_to_clade[gca]
        bgcome_size = float(ls[4])/1000000.0
        comp_bgc_count = float(ls[3])
        strict_nrps_or_pksome_size = float(ls[-1])/1000000.0
        
        if clade == 'clade_1':
            clade_1_bgcomes.append(bgcome_size*1000000.0)
        else:
            clade_2_bgcomes.append(bgcome_size*1000000.0)

        print('\t'.join([gca, clade, str(comp_bgc_count), str(bgcome_size),  str(strict_nrps_or_pksome_size), grid]))

sys.stderr.write(str(stats.ranksums(clade_1_bgcomes, clade_2_bgcomes, alternative='greater')[1]) + '\n')
