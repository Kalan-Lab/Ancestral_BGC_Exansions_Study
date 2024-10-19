import os
import sys

print('DATASET_MULTIBAR')
print('SEPARATOR TAB')
print('DATASET_LABEL\tBGC_Size_Relaxed_and_Strict')
print('COLOR\t#000000')
print('FIELD_LABELS\tStrict\tRelaxed')
print('FIELD_COLORS\t#747474\t#b5b5b3')
print('DATA')

with open('../AntiSMASH_Stats.txt') as oaf:
    for i, line in enumerate(oaf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        print(ls[1] + '\t' + ls[-1] + '\t' + str(int(ls[4])-int(ls[-1])))

