import os
import sys
from collections import defaultdict

select_genomes_file = 'Select_Leaves.txt'
cluster_file = 'zol_dom_Results/Determine_Orthogroups/Orthogroups.tsv'

select_genomes = set([])
with open(select_genomes_file) as osgf:
    for i, line in enumerate(osgf):
        line = line.strip()
        select_genomes.add(line)

strep = set(['BGC0001785', 'GCA_024611995', 'GCA_001746425', 'GCA_004364395', 'GCA_004794975', 'GCA_017948485', 'GCA_001514035', 'GCA_022696125', 'GCA_020907265', 'GCA_014648855', 'GCA_001660045', 'GCA_900105265', 'BGC0000106', 'GCA 001013905', 'BGC0001287'])
amyco = set(['GCA_017308975', 'GCA_000716785', 'BGC0002009', 'BGC0001759', 'GCA_013364075', 'GCA_900105945', 'GCA_000695625', 'BGC0000136', 'GCA_001742805', 'GCA_000454025', 'GCA_000700945', 'GCA_000196835', 'GCA_000282715', 'GCA_000696405', 'GCA_000220945'])
micro = set(['GCA_020884795', 'GCA_900090265', 'GCA_011765705', 'GCA_011765735', 'BGC0000137', 'GCA_000384275', 'GCA_011765725'])

genomes = []
genome_ogs = defaultdict(set)
with open(cluster_file) as ocf:
    for i, line in enumerate(ocf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: 
            genomes = [x.strip().split('.')[0] for x in ls[1:]]
        else:
            og = ls[0]
            main_og = False
            ahba_og = False
            genomes_with = set([])
            for j, lts in enumerate(ls[1:]):
                g = genomes[j]
                if g in select_genomes:
                    for lt in lts.split(', '):
                        lt = lt.strip()
                        if lt == '': continue
                        genome_ogs[g].add(og)
           
taxa = [strep, amyco, micro]
taxa_name = ['strep', 'amyco', 'micro']

print('comparison\tjaccard_index')
for i, t1 in enumerate(taxa):
    tn1 = taxa_name[i]
    for j, t2 in enumerate(taxa):
        tn2 = taxa_name[j]
        if i > j: continue
        for k, g1 in enumerate(t1):
            g1_ogs = genome_ogs[g1]
            for l, g2 in enumerate(t2):
                if k >= l: continue
                g2_ogs = genome_ogs[g2]
                g1_g2_jci = len(g1_ogs.intersection(g2_ogs))/len(g1_ogs.union(g2_ogs))
                print(tn1 + ' vs. ' + tn2 + '\t' + str(g1_g2_jci))
