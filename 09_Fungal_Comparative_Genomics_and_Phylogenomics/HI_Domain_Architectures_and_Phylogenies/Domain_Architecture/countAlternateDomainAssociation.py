import os
import sys
from collections import defaultdict

alt_domains = defaultdict(set)
focal_dom = sys.argv[2]

with open(sys.argv[1]) as oad:
    for line in oad:
        line = line.strip()
        ls = line.split('\t')
        doms = line.split('\t')[1].split('|')
        if len(doms) == 1:
            alt_domains['standalone'].add(ls[0])
        else:
            for dom in doms:
                if dom != focal_dom:
                    alt_domains[dom].add(ls[0])

for dom in alt_domains:
    print(dom + '\t' + str(len(alt_domains[dom])))
