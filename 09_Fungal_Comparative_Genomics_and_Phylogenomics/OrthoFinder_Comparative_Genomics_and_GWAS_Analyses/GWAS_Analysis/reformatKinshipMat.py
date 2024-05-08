import os
import sys

with open(sys.argv[1]) as opf:
    for i, line in enumerate(opf):
        line = line.strip('\n')
        ls = line.split('\t')
        printlist = []
        if i == 0:
            printlist.append(ls[0])
            for val in ls[1:]:
                printlist.append('_'.join(val.split('_')[-2:]))
        else:
            printlist.append('_'.join(ls[0].split('_')[-2:]))
            printlist += ls[1:]
        print('\t'.join(printlist))
