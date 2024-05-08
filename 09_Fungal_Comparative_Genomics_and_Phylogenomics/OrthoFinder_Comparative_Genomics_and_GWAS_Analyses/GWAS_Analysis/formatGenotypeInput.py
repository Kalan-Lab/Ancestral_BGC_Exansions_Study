import os
import sys

orthofinder_tsv = sys.argv[1]

with open(orthofinder_tsv) as oot:
    for i, line in enumerate(oot):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: 
            print(line)
        else:
            printlist = [ls[0]]
            for val in ls[3:]:
                if val.strip() != '':
                    printlist.append('1')
                else:
                    printlist.append('0')
            print('\t'.join(printlist))
