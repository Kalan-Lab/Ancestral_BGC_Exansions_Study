import os
import sys

as_dir = os.path.abspath(sys.argv[1]) + '/'

for f in os.listdir(as_dir):
    if f.endswith('.gbk') and '.region' in f:
        scaff = f.split('.region')[0]
        start = None
        end = None
        with open(as_dir + f) as oaf:
            for line in oaf:
                line = line.strip()
                if line.startswith('Orig. start  ::'): start = line.split()[-1]
                elif line.startswith('Orig. end    ::'): end = line.split()[-1]
        print(f.split('.gbk')[0] + '\t' + scaff + '\t' + start + '\t' + end)
