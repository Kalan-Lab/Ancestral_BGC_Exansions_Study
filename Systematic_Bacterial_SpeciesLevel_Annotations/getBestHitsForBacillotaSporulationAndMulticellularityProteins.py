import os
import sys
from collections import defaultdict

cog_bs_file = 'COG_Bitscore_Cutoffs.txt'
results_file = 'Diamond_Results.txt'

cog_bitscore_cutoffs = {}
with open(cog_bs_file) as ocbf:
    for line in ocbf:
        line = line.strip()
        ls = line.split('\t')
        cog_bitscore_cutoffs[ls[0]] = float(ls[1])

top_hits = defaultdict(lambda: [[], 0.0])
with open(results_file) as orf:
    for line in orf:
        line = line.strip()
        qseqid, sseqid, pident, length, evalue, bitscore, qcovhsp, scovhsp = line.split('\t')
        bitscore = float(bitscore)
        cog = qseqid.split('|')[0]
        if bitscore > cog_bitscore_cutoffs[cog]:
            if bitscore > top_hits[sseqid][1]:
                top_hits[sseqid] = [[cog], bitscore]
            elif bitscore == top_hits[sseqid][1]:
                top_hits[sseqid][0].append(cog)

for prot in top_hits:
    print(prot + '\t' + '; '.join(top_hits[prot][0]) + '\t' + str(top_hits[prot][1]))
