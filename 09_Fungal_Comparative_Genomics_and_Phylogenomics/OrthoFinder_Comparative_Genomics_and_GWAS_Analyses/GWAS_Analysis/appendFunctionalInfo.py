import os
import sys
from collections import defaultdict
from operator import itemgetter

sighog_file = 'Significant_Associations_Dikarya.txt'
eggnog_file = 'out.emapper.annotations'

cog_cats = defaultdict(lambda: defaultdict(int))
desc_cats = defaultdict(lambda: defaultdict(int))
with open(eggnog_file) as oef:
    for line in oef:
        if line.startswith('#'): continue
        line = line.strip()
        ls = line.split('\t')
        hog = ls[0].split('|')[0]
        cogs = ls[6]
        desc = ls[7]
        for cog in cogs:
            cog_cats[hog][cog] += 1
        desc_cats[hog][desc] += 1

for i, line in enumerate(open(sighog_file)):
    line = line.strip()
    if i == 0: 
        print(line)
    else:
        ls = line.split('\t')
        hog = ls[0].split('N0.')[1]
        
        desc_string = []
        for tup in sorted(desc_cats[hog].items(), key=itemgetter(1), reverse=True):
            desc_string.append(tup[0] + '=' + str(tup[1]))

        cog_string = []
        for tup in sorted(cog_cats[hog].items(), key=itemgetter(1), reverse=True):
            cog_string.append(tup[0] + '=' + str(tup[1]))
        print('\t'.join([hog, '; '.join(desc_string), '; '.join(cog_string)]))
