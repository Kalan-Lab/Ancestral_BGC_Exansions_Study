import os
import sys

hog_file = sys.argv[1]
for line in open(hog_file):
    line = line.strip('\n')
    ls = line.split('\t')
    print(ls[0] + '\t' + '\t'.join(ls[3:]))

