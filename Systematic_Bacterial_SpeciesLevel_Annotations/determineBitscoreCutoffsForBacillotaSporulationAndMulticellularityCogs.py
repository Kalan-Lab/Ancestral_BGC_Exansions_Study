import os
import sys
from collections import defaultdict

cog_min_bitscores = defaultdict(lambda: 1.0e100)
with open('Reflexive_Diamond_Results.txt') as ordrf:
    for line in ordrf:
        line = line.strip()
        ls = line.split('\t')
        cog1 = ls[0].split('|')[0]
        cog2 = ls[1].split('|')[0]
        if cog1 == cog2:
            bitscore = float(ls[3])
            if bitscore < cog_min_bitscores[cog1]:
                cog_min_bitscores[cog1] = bitscore

for cog in cog_min_bitscores:
    print(cog + '\t' + str(cog_min_bitscores[cog]))

