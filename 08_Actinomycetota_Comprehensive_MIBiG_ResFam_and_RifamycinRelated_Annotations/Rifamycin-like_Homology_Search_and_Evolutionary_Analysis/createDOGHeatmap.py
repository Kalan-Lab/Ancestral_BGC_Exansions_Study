import os
import sys
from collections import defaultdict


main_rif_genes = set(['AABH_000021', 'AABH_000022', 'AABH_000023', 'AABH_000024', 'AABH_000025']) # rifABCDE
ahba_genes = set(['AABH_000028', 'AABH_000029', 'AABH_000030', 'AABH_000031', 'AABH_000032', 'AABH_000033', 'AABH_000034', 'AABH_000060']) # rifGHIKLMN-J

select_genomes_file = 'Select_Leaves.txt'
cluster_file = 'zol_dom_Results/Determine_Orthogroups/Orthogroups.tsv'

select_genomes = set([])
with open(select_genomes_file) as osgf:
    for i, line in enumerate(osgf):
        line = line.strip()
        select_genomes.add(line)

genomes = []
og_genomes = defaultdict(set)
main_ogs = set([])
ahba_ogs = set([])
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
                        genomes_with.add(g)
                        gene = lt.split('|')[1]
                        if gene in main_rif_genes:
                            main_og = True
                        elif gene in ahba_genes:
                            ahba_og = True
            
            if len(genomes_with) < 2: continue
            for g in genomes_with:
                og_genomes[og].add(g)

            if main_og == True and ahba_og == True:
                print('ISSUE')
                print(line)
                sys.exit(1)

            if main_og:
                main_ogs.add(og)
            elif ahba_og:
                ahba_ogs.add(og)

print('DATASET_HEATMAP')
print('SEPARATOR TAB')
print('DATASET_LABEL\tDOGs')
print('COLOR\t#000000')
print('FIELD_LABELS\t' + '\t'.join(sorted(og_genomes)))
print('DATA')
for g in select_genomes:
    row = [g]
    for og in sorted(og_genomes):
        if g in og_genomes[og]:
            if og in main_ogs:
                row.append('3')
            elif og in ahba_ogs:
                row.append('2')
            else:
                row.append('1')
        else:
            row.append('NA')
    print('\t'.join(row))
