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

fam_gca_bgcs = defaultdict(list)
fam_gca_bgcs_counts = defaultdict(lambda: defaultdict(int))
for gca in gca_bgc_prots:
    for bgc in gca_bgc_prots[gca]:
        gca_bgc_key_prot_found = len(gca_bgc_prots)
        gca_bgc_key_prot_found_prop = gca_bgc_key_prot_found/float(bgc_keyprot_counts[bgc])
        if gca_bgc_key_prot_found_prop >= 0.5:
            fam = gca_to_fam[gca]
            fam_gca_bgcs[fam + '|' + gca].append(bgc)
            fam_gca_bgcs_counts[fam + '|' + gca][bgc] += 1

for fg in fam_gca_bgcs:
    f, g = fg.split('|')
    print(f + '\t' + g + '\t' + str(', '.join(sorted(fam_gca_bgcs[fg]))))
