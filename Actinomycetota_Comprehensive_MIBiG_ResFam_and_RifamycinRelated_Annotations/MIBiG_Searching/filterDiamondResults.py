import os
import sys

diamond_result_file = "Search_Full_Actinos_for_MIBiG_Key_Prots.txt" # Not included because it is 12GB gzipped!

with open(diamond_result_file) as odrf:
    for i, line in enumerate(odrf):
        if (i % 100000) == 0: sys.stderr.write('processed %s\n' % i) 
        line = line.strip()
        query, hit, pid, length, evalue, bitscore, qcovhsp, scovhsp = line.split('\t')
        pid = float(pid)
        qcovhsp = float(qcovhsp)
        scovhsp = float(scovhsp)
        if pid >= 95.0 and qcovhsp >= 70.0 and scovhsp >= 70.0:
            print(line)
