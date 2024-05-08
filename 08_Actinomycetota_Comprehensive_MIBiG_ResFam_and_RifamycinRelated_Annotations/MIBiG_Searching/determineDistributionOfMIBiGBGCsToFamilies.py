import os
import sys
from collections import defaultdict
from Bio import SeqIO

gca_to_tax_file = 'GCA_to_Taxonomy.txt'
mibig_results_file = 'Search_Full_Actinos_for_MIBiG_Key_Prots.Filtered.txt'
tree_fam_file = 'FamilyTree_Represented.txt'
bgc_fasta_file = 'MIBiG_KeyBiosynthetic_Proteins.faa'

bgc_keyprot_counts = defaultdict(int)
with open(bgc_fasta_file) as obff:
    for rec in SeqIO.parse(obff, 'fasta'):
        bgc_keyprot_counts[rec.id.split('|')[0]] += 1

gca_to_fam = {}
fam_to_cla = {}
fam_gca_counts = defaultdict(int)
with open(gca_to_tax_file) as ogttf:
    for line in ogttf:
        line = line.strip()
        ls = line.split('\t')
        fam = ls[1].split('f__')[1].split(';g__')[0]
        cla = ls[1].split('c__')[1].split(';o__')[0]
        gca_to_fam[ls[0]] = fam
        fam_to_cla[fam] = cla
        fam_gca_counts[fam] += 1
    
gca_bgc_prots = defaultdict(lambda: defaultdict(set))
fam_bgc_ids = defaultdict(set)
with open(mibig_results_file) as omrf:
    for line in omrf:
        line = line.strip()
        ls = line.split('\t')
        bgc = ls[0].split('|')[0]
        gca = '_'.join(ls[1].split('_')[:2]).split('.')[0]
        gca_bgc_prots[gca][bgc].add(ls[0])

for gca in gca_bgc_prots:
    for bgc in gca_bgc_prots[gca]:
        gca_bgc_key_prot_found = len(gca_bgc_prots[gca][bgc])
        gca_bgc_key_prot_found_prop = gca_bgc_key_prot_found/float(bgc_keyprot_counts[bgc])
        if gca_bgc_key_prot_found_prop >= 0.5:
            fam = gca_to_fam[gca]
            fam_bgc_ids[fam].add(bgc)

fam_accounted = set([])
with open(tree_fam_file) as otff:
    for line in otff:
        fam = line.strip()
        print(fam + '\t' + fam_to_cla[fam] + '\tIn Phylogeny\t' +  str(len(fam_bgc_ids[fam])) + '\t' + str(fam_gca_counts[fam]))
        fam_accounted.add(fam)

for fam in fam_bgc_ids:
    if not fam in fam_accounted:
        print(fam + '\t' + fam_to_cla[fam] + '\tNot In Phylogeny\t' + str(len(fam_bgc_ids[fam])) + '\t' + str(fam_gca_counts[fam]))
