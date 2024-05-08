import os
import sys
from collections import defaultdict

rae_distribution_file = 'GCA_RAE_Palindromes_Listing.txt'
rif_distribution_file = sys.argv[1] 

gcas_with_rae = set([])
with open(rae_distribution_file) as ordf:
    for line in ordf:
        line = line.strip()
        ls = line.split('\t')
        gcas_with_rae.add(ls[0])

with open(rif_distribution_file) as ordf:
    for line in ordf:
        line = line.strip()
        gca = line.split('.')[0]
        rae_found = 'False'
        if gca in gcas_with_rae:
            rae_found = 'True'
        print(line + '\t' + rae_found)
