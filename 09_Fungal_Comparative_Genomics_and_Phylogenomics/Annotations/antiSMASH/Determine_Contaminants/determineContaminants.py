import os
import sys
from collections import defaultdict

query_hits = defaultdict(int)
with open('DIAMOND_blastp_to_bacNR_Results.txt') as orf:
    for line in orf:
        line = line.strip()
        ls = line.split('\t')
        pid = float(ls[2])
        qcov = float(ls[-1])
        alglen = float(ls[3])
        evalue = float(ls[-3])
        if pid >= 80.0 and qcov >= 80.0: 
            query_hits[ls[0]] += 1
            
for query in query_hits:
    print(query.split('|')[0] + '\t' + query.split('|')[2] + '\t' + str(query_hits[query]))
