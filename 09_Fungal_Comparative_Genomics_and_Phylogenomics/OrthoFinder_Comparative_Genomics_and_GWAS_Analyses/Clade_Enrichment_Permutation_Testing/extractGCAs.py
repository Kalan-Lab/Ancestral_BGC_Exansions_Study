import os
import sys

with open(sys.argv[1]) as of:
    for line in of:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        print(gca)
