import os
import sys
from collections import defaultdict

rae_hit_dir = 'BLASTn_Results/'
gca_to_tax_file = 'GCA_to_Taxonomy.txt'
tree_fam_file = 'FamilyTree_Represented.txt'

gca_to_fam = {}
fam_to_cla = {}
with open(gca_to_tax_file) as ogttf:
    for line in ogttf:
        line = line.strip()
        ls = line.split('\t')
        fam = ls[1].split('f__')[1].split(';g__')[0]
        cla = ls[1].split('c__')[1].split(';o__')[0]
        gca_to_fam[ls[0]] = fam
        fam_to_cla[fam] = cla
    
fam_hits = defaultdict(int)
fam_total = defaultdict(int)
for f in os.listdir(rae_hit_dir):
    rae_hit_file = rae_hit_dir + f
    hit = 0
    with open(rae_hit_file) as orhf:
        for line in orhf:
            line = line.strip()
            ls = line.split('\t')
            if float(ls[2]) == 100.0 and float(ls[-1]) == 100.0:
                hit += 1
    gca = '_'.join(f.split('_')[:2]).split('.')[0]
    fam = gca_to_fam[gca]
    fam_total[fam] += 1
    if hit >= 1:
        fam_hits[fam] += 1

clade1 = set([])
clade2 = set([])
with open('Clade_1_Families.txt') as oc1f:
    for line in oc1f:
        line = line.strip()
        clade1.add(line)

with open('Clade_2_Families.txt') as oc2f:
    for line in oc2f:
        line = line.strip()
        clade2.add(line)

print('Family\tClassGroup\tClass\tInPhylo\tWithOrWithout\tProportion\tProportionWith\tTotalGenomes')
fam_accounted = set([])
with open(tree_fam_file) as otff:
    for line in otff:
        fam = line.strip()
        fam_group = '2Actinomycetia'
        if fam_to_cla[fam] == 'Actinomycetia':
            if fam in clade1:
                fam_group = '3Actinomycetia - Clade-1'
            elif fam in clade2:
                fam_group = '4Actinomycetia - Clade-2'
        else:
            fam_group = '5Other Actinomycetota'
        print(fam + '\t' + fam_group + '\t' + fam_to_cla[fam] + '\tIn_Phylogeny\t7With\t' + str(fam_hits[fam]/float(fam_total[fam])) + '\t' + str(fam_hits[fam]/float(fam_total[fam])) + '\t' + str(fam_total[fam])) 
        print(fam + '\t' + fam_group + '\t' + fam_to_cla[fam] + '\tIn_Phylogeny\t6Without\t' + str(1.0-(fam_hits[fam]/float(fam_total[fam]))) + '\t' + str(fam_hits[fam]/float(fam_total[fam])) + '\t' + str(fam_total[fam]))
        print(fam + '\t' + fam_group + '\t' + fam_to_cla[fam] + '\tNot_In_Phylogeny\t' + fam_group + '\t0.1\t' + str(fam_hits[fam]/float(fam_total[fam])) + '\t' + str(fam_total[fam]))
        fam_accounted.add(fam)

for fam in fam_total:
    if not fam in fam_accounted:
        fam_group = '2Actinomycetia'
        if fam_to_cla[fam] == 'Actinomycetia':
            if fam in clade1:
                fam_group = '3Actinomycetia - Clade-1'
            elif fam in clade2:
                fam_group = '4Actinomycetia - Clade-2'
        else:
            fam_group = '5Other Actinomycetota'
        print(fam + '\t' + fam_group + '\t' + fam_to_cla[fam] + '\tNot_In_Phylogeny\t7With\t' + str(fam_hits[fam]/float(fam_total[fam])) + '\t' + str(fam_hits[fam]/float(fam_total[fam])) + '\t' + str(fam_total[fam]))
        print(fam + '\t' + fam_group + '\t' + fam_to_cla[fam] + '\tNot_In_Phylogeny\t6Without\t' + str(1.0-(fam_hits[fam]/float(fam_total[fam]))) + '\t' + str(fam_hits[fam]/float(fam_total[fam])) + '\t' + str(fam_total[fam]))
        print(fam + '\t' + fam_group + '\t' + fam_to_cla[fam] + '\tNot_In_Phylogeny\t' + fam_group + '\t0.1\t' + str(fam_hits[fam]/float(fam_total[fam])) + '\t' + str(fam_total[fam]))
