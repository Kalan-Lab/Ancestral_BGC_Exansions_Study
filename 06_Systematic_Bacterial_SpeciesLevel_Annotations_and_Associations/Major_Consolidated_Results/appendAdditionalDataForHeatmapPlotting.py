import os
import sys
from collections import defaultdict

raw_data_file = sys.argv[1]

factor_order = defaultdict(int)
with open(raw_data_file) as ordf:
    for line in ordf:
        line = line.strip()
        ls = line.split('\t')
        factor_order[ls[1]] += float(ls[2])

phylum_order = {'General': '6', 'Actinomycetota': '5', 'Bacillota': '2', 'Cyanobacteriota': '3', 'Myxococcota': '4', 'Pseudomonadota': '1'}
print('Phylum_Order\tPhylum\tFactor\tFactor_Order\tSpearman_Correlation\tPval\tPval_Category')
with open(raw_data_file) as ordf:
    for i, line in enumerate(ordf):
        line = line.strip()
        ls = line.split('\t')
        pval = float(ls[-1])
        pval_cat = '1'
        if pval < 1e-3:
            pval_cat = '2'
        if pval < 1e-5:
            pval_cat = '3'
        if pval < 1e-10:
            pval_cat = '4'
        phylum_ord = phylum_order[ls[0]]
        print('\t'.join([phylum_ord, ls[0], ls[1], str(factor_order[ls[1]]), ls[2], ls[-1], pval_cat]))
