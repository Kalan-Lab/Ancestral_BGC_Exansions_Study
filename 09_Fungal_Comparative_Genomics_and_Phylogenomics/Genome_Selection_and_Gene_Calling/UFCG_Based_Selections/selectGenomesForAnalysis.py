import os
import sys
from collections import defaultdict

gca_with_genomes = set([])
genomes_dir = 'Genomes/'
for f in os.listdir(genomes_dir):
    gca = '_'.join(f.split('_')[:2])
    gca_with_genomes.add(gca)

abyss_file = 'Abyss_Fac_Stats.txt'
gca_n50s = {}
with open(abyss_file) as oaf:
    for i, line in enumerate(oaf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca = '_'.join(ls[-1].split('/')[-1].split('_')[:2])
        n50 = ls[5]
        gca_n50s[gca] = n50

family_listing_file = 'NonAscomycota_Classes.txt' #'Families_of_Interest.txt'
select_families = set([])
with open(family_listing_file) as oflf:
    for line in oflf:
        line = line.strip()
        select_families.add(line)

ufcg_dir = 'Species_UFCG/'
family_reps = defaultdict(list)
for f in os.listdir(ufcg_dir):
    gca = '_'.join(f.split('_')[:2])
    with open(ufcg_dir + f) as ouf:
        for line in ouf:
            line = line.strip()
            if line.startswith('"taxonomy":'):
                fam = line.split(';')[2] # [4]
                if fam in select_families and gca in gca_with_genomes:
                    n50 = gca_n50s[gca]
                    family_reps[fam].append([gca, n50])

for fam in family_reps:
    max_n50 = max([x[1] for x in family_reps[fam]])
    for gca_info in family_reps[fam]:
        gca, n50 = gca_info
        if n50 == max_n50:
            print(fam + '\t' + gca)
