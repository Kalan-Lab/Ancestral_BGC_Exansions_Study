import os
import sys
from collections import defaultdict
import json
import math

"""
COLUMN INFORMATION:
----------------------------------------
0       npaid
1       compound_id
2       compound_names
3       compound_molecular_formula
4       compound_molecular_weight
5       compound_accurate_mass
6       compound_m_plus_h
7       compound_m_plus_na
8       compound_inchi
9       compound_inchikey
10      compound_smiles
11      compound_cluster_id
12      compound_node_id
13      origin_type
14      genus
15      origin_species
16      original_reference_author_list
17      original_reference_year
18      original_reference_issue
19      original_reference_volume
20      original_reference_pages
21      original_reference_doi
22      original_reference_pmid
23      original_reference_title
24      original_reference_type
25      original_journal_title
26      reassignment_dois
27      synthesis_dois
28      mibig_ids
29      gnps_ids
30      npmrd_id
31      npatlas_url
"""

# This script is a little more involved - because we were attempting to additionally gather 
# annotation (npclassifier superpathway information) but there are too many types and we 
# resorted to just plotting the counts 

npa_json = 'NPAtlas_download_2023_06.json' # The JSON file with in depth info on metabolites from NPAtlas 
npa_tsv = 'NPAtlas_download_2023_06.tsv' # The TSV file describing major info for metabolites from NPAtlas
conc_phy = 'Concatenated_Alignment.Renamed.phy' # Just the Phylip formatted alignment used during MCMCTree construction of the timetree (to get the species names)
gtdb_col = 'GTDB_Taxa_Column.txt' # the 16th column corresponding to GTDB R214 taxonomies of genomes in the release.

np_to_pathway = {}
with open(npa_json) as onpaj:
    npaj = json.load(onpaj)
    for comp in npaj:
        comp_id = comp['npaid']
        npclassifier_superpathways = []
        try:
            npclassifier_superpathways = comp['npclassifier']['superclass_results']
        except:
            pass
        if len(npclassifier_superpathways) == 0 or len(npclassifier_superpathways) > 1:
            npsp = 'Multiple/Other'
        else:
            npsp = npclassifier_superpathways[0]
        np_to_pathway[comp_id] = npsp 

orders_of_interest = set([])
with open(conc_phy) as ocf:
    for i, line in enumerate(ocf):
        if i == 0: continue
        line = line.strip()
        ls = line.split()
        orders_of_interest.add(ls[0])

genus_to_order = {}
with open(gtdb_col) as ogtcf:
    for i, line in enumerate(ogtcf):
        if i == 0: continue 
        line = line.strip()
        genus = line.split(';g__')[1].split(';s__')[0]
        order = line.split(';o__')[1].split(';f__')[0]
        genus_to_order[genus] = order

order_nodes = defaultdict(lambda: defaultdict(set))
with open(npa_tsv) as of:
    for i, line in enumerate(of):
        if i == 0: continue
        line = line.strip('\n')
        ls = line.split('\t')
        npaid = ls[0]
        node = ls[12]
        genus = ls[14]
        if not genus in genus_to_order: continue
    order = genus_to_order[genus]
        order_nodes[order][node].add(np_to_pathway[npaid])

order_node_type_counts = defaultdict(lambda: defaultdict(int))
for order in orders_of_interest:
    for node in order_nodes[order]:
        if len(order_nodes[order][node]) > 1:
            order_node_type_counts[order]['Multiple/Other'] += 1
        else:
            order_node_type_counts[order][list(order_nodes[order][node])[0]] += 1


"""
colors = []
with open('Colors_for_Major_Superpathways.txt') as ocf:
    for line in ocf:
        line = line.strip()
        colors.append(line)

annots_in_order = []
with open('Major_NPClassifier_SuperPathways.txt') as omf:
    for line in omf:
        line = line.strip()
        annots_in_order.append(line)
"""

print('DATASET_SIMPLEBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tNPAtlas Distinct Nodes (Simple)')
print('DATA')

for order in orders_of_interest:
    count = 0
    for nt in order_node_type_counts[order]:
        count += order_node_type_counts[order][nt]
    if count == 0: continue
    print(order + '\t' + str(math.log(count, 10)))

"""
print('DATASET_MULTIBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tMajor Annotations')
print('FIELD_LABELS\tMultiple/Other\t' + '\t'.join(annots_in_order))
print('FIELD_COLORS\t#4f4f4f\t' + '\t'.join(colors))
print('DATA')

for order in orders_of_interest:
    other_count = 0
    printlist = []
    for nt in order_node_type_counts[order]:
        if not nt in annots_in_order:
            other_count += order_node_type_counts[order][nt]
    for nt in annots_in_order:
        printlist.append(order_node_type_counts[order][nt])
    printlist = [order, other_count] + printlist
    print('\t'.join([str(x) for x in printlist]))
"""
