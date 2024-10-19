import os
import sys

for line in open('uniprot_blastp_results.txt'):
    line = line.strip()
    ls = line.split('\t')
    if float(ls[2]) >= 90.0:
        print(line)
