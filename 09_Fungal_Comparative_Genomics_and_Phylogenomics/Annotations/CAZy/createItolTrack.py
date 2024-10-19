import os
import sys
from collections import defaultdict

gca_to_name = {}
with open('../../Overview_File.txt') as ofgg:
    for line in ofgg:
        line = line.strip()
        ls = line.split('\t')
        gca_to_name[ls[0]] = ls[1]

print('DATASET_SIMPLEBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tCAZy_Homologs_Distinct')
print('COLOR\t#000000')
print('DATA')

gca_counts = defaultdict(set)
with open('CAZy_Summary.txt') as othb:
    for line in othb:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        print(gca_to_name[gca] + '\t' + str(ls[1]))
