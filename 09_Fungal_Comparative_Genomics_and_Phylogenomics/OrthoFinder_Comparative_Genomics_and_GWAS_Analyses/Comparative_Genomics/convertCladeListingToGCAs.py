import os
import sys

reps = set([x.strip().split('\t')[0] for x in open('Representatives.txt').readlines()])

with open(sys.argv[1]) as of:
    for line in of:
        line = line.strip()
        gca = '_'.join(line.split('_')[1:]).split('.')[0].replace('GCF_', 'GCA')
        if gca in reps:
            print(gca)
