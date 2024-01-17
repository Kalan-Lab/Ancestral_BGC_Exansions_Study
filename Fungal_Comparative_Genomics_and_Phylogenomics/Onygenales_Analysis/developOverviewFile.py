import os
import sys
from collections import defaultdict

gca_to_taxid_file = 'Onygenales_GCAs.txt'
best_annots_file = 'Best_Available_Annotations.txt'

gca_to_taxonomy = {}
with open(gca_to_taxid_file) as ogttf:
    for line in ogttf:
        line = line.strip()
        ls = line.split('\t')
        gca_to_taxonomy[ls[1].split('.')[0]] = ls[0]

print('gca\ttaxonomy\tannotation_method\tversion_info\tncbi_annot_avail\tfunannotate_annot_avail\tfasta_file\tgff_file')
with open(best_annots_file) as obaf:
    for line in obaf:
        line = line.strip()
        ls = line.split('\t')
        print('\t'.join([ls[0], gca_to_taxonomy[ls[0]]] + ls[1:]))
