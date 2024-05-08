import os
import sys
from collections import defaultdict
from operator import itemgetter

sighog_file = 'Significant_HOGs.txt'
eggnog_file = 'out.emapper.annotations.sorted'
bgc_file = 'All_BGC_Proteins.with_Labels.txt'

cog_cats = defaultdict(lambda: defaultdict(int))
desc_cats = defaultdict(lambda: defaultdict(int))
with open(eggnog_file) as oef:
    for line in oef:
        if line.startswith('#'): continue
        line = line.strip()
        ls = line.split('\t')
        hog = ls[0].split('N0.')[1].split('|')[0]
        cogs = ls[6]
        desc = ls[7]
        for cog in cogs:
            cog_cats[hog][cog] += 1
        desc_cats[hog][desc] += 1

bgc_region_found = set([])
key_bgc_hog = set([])
with open(bgc_file) as obf:
    for line in obf:
        line = line.strip()
        ls = line.split('\t')
        if not 'N0' in ls[0]: continue
        hog = ls[0].split('N0.')[1]
        bgc_region_found.add(hog)
        if ls[-1] == 'True':
            key_bgc_hog.add(hog)

for i, line in enumerate(open(sighog_file)):
    line = line.strip()
    if i == 0: 
        print(line)
    else:
        ls = line.split('\t')
        hog = ls[1].split('N0.')[1]
        
        found_in_bgc = False
        cat = 'NA'
        if hog in bgc_region_found:
            found_in_bgc = True
            if hog in key_bgc_hog:
                cat = 'Key BGC protein'
        desc_string = []
        for tup in sorted(desc_cats[hog].items(), key=itemgetter(1), reverse=True):
            desc_string.append(tup[0] + '=' + str(tup[1]))

        cog_string = []
        for tup in sorted(cog_cats[hog].items(), key=itemgetter(1), reverse=True):
            cog_string.append(tup[0] + '=' + str(tup[1]))
        print('\t'.join([hog, '; '.join(desc_string), '; '.join(cog_string), cat, str(found_in_bgc), ls[0]] + ls[2:]))
