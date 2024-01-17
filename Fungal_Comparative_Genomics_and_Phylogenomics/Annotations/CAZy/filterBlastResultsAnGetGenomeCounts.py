import os
import sys
from collections import defaultdict

top_hits = defaultdict(lambda: [[], 0.0]) # [0.0], [0.0]])
with open('CAZy_Blastp_Results.txt') as obr:
    for line in obr:
        line = line.strip()
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qcovhsp, scovhsp = line.split('\t')
        pident = float(pident)
        qcovhsp = float(qcovhsp)
        scovhsp = float(scovhsp)
        bitscore = float(bitscore)
        if pident < 30.0 or qcovhsp < 70.0 or scovhsp < 70.0: continue
        if top_hits[qseqid][1] < bitscore:
            top_hits[qseqid] = [[sseqid], bitscore]# pident, qcovhsp]
        elif top_hits[qseqid][1] == bitscore:
            top_hits[qseqid][0].append(sseqid)

for s in top_hits:
    print(s + '\t' + str(top_hits[s][1]) + '\t' + '\t'.join(top_hits[s][0]))

