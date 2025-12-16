import os
import sys

print('DATASET_COLORSTRIP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tConsidered')
print('DATA')
with open('RifK_to_FaiGC.txt') as orf:
    for line in orf:
        line = line.strip()
        ls = line.split("\t")
        print(ls[0] + '\t#000000\tconsidered')
