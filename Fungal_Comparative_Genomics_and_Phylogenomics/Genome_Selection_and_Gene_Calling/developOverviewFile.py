import os
import sys
from collections import defaultdict

gca_to_taxid_file = 'GCA_to_Taxid.txt'
taxid_to_name_file = 'Taxids_to_Taxonomy_Info.txt'
best_annots_file = 'Best_Available_Annotations.txt'

taxid_to_taxonomy = {}
with open(taxid_to_name_file) as otnf:
    for line in otnf:
        line = line.strip()
        ls = line.split('\t')
        taxid_to_taxonomy[ls[0]] = ls[1]

gca_to_taxonomy = {}
with open(gca_to_taxid_file) as ogttf:
    for i, line in enumerate(ogttf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca_to_taxonomy[ls[0]] = taxid_to_taxonomy[ls[2]]

print('gca\ttaxonomy\tannotation_method\tversion_info\tncbi_annot_avail\tfunannotate_annot_avail\tfasta_file\tgff_file')
with open(best_annots_file) as obaf:
    for line in obaf:
        line = line.strip()
        ls = line.split('\t')
        print('\t'.join([ls[0], gca_to_taxonomy[ls[0]]] + ls[1:]))
