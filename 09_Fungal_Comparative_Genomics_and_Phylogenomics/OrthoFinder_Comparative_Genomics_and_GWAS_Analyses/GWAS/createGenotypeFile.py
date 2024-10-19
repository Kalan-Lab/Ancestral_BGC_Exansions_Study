import os
import sys
from collections import defaultdict

orthogroup_count_tsv_file = 'Orthogroups.tsv'
with open(orthogroup_count_tsv_file) as octf:
    for i, line in enumerate(octf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            print('\t'.join(ls))
        else:
            printrow = [ls[0]]
            for val in ls[1:]:
                if val == '':
                    printrow.append('0')
                else:
                    printrow.append('1')
            print('\t'.join(printrow))
