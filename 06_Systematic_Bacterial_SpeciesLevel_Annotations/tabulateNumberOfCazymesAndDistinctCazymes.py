import os
import sys
from collections import defaultdict

# thresholds from : https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/readme.txt

evalue_threshold = 1e-18 
coverage_threshold = 0.35

genome_caz = defaultdict(set)
genome_caz_total = defaultdict(set)
with open('PyHMMER_Results.txt') as oprf:
    for line in oprf:
        line = line.strip()
        ls = line.split('\t')
        cov = float(ls[6])
        evalue = float(ls[4])
        if evalue < evalue_threshold and cov > coverage_threshold:
            caz = ls[0].split('.hmm')[0]
            genome = ls[2].split('|')[0]
            genome_caz[genome].add(caz)
            genome_caz_total[genome].add(ls[2])

for genome in genome_caz_total:
    print(genome + '\t' + str(len(genome_caz_total[genome])) + '\t' + str(len(genome_caz[genome])))
