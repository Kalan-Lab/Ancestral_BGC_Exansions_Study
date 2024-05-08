import os
import sys
from collections import defaultdict

gca_to_tax_file = 'GCA_to_Taxonomy.txt'
mibig_results_file = 'Search_Full_Actinos_for_MIBiG_Key_Prots.Filtered.txt'
clustered_mibig_file = 'Slclust_Clustered.with_Singletons.txt'
tree_fam_file = 'FamilyTree_Represented.txt'

prot_to_id = {}
with open(clustered_mibig_file) as ocmf:
    for i, line in enumerate(ocmf):
        line = line.strip()
        ls = line.split('\t')
        for prot in ls:
            prot_to_id[prot] = i

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
    
fam_prot_ids = defaultdict(set)
with open(mibig_results_file) as omrf:
    for line in omrf:
        line = line.strip()
        ls = line.split('\t')
        prot_id = prot_to_id[ls[0]] 
        gca = '_'.join(ls[1].split('_')[:2]).split('.')[0]
        fam = gca_to_fam[gca]
        fam_prot_ids[fam].add(prot_id)

fam_accounted = set([])
with open(tree_fam_file) as otff:
    for line in otff:
        fam = line.strip()
        print(fam + '\t' + fam_to_cla[fam] + '\tIn Phylogeny\t' + str(len(fam_prot_ids[fam])))
        fam_accounted.add(fam)

for fam in fam_prot_ids:
    if not fam in fam_accounted:
        print(fam + '\t' + fam_to_cla[fam] + '\tNot In Phylogeny\t' + str(len(fam_prot_ids[fam])))
