import os
import sys

gca_to_id = {}
with open('AntiSMASH_Stats.Updated_WithoutContaminants.txt') as occf:
    for i, line in enumerate(occf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca_to_id[ls[0]] = ls[1]

print('DATASET_SIMPLEBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tCAZy_Total')
print('COLOR\t#000000')
print('DATA')
with open('CAZy_Counts_HMM_Based.txt') as of:
    for line in of:
        line = line.strip()
        ls = line.split('\t')
        total = int(ls[1])
        print(gca_to_id[ls[0]] + '\t' + str(total))
