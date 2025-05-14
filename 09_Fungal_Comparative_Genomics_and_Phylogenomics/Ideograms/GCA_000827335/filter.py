import os
import sys

keep = set([])
with open('BGC_Lengths.txt') as of:
    for line in of:
        line = line.strip()
        ls = line.split()
        if int(ls[2]) >= 1000000:
            keep.add(ls[1])

with open(sys.argv[1]) as of:
    for i, line in enumerate(of):
        line = line.strip()
        if i == 0: 
            print(line)
            continue
        ls = line.split()
        if ls[0] in keep:
            print(line)

