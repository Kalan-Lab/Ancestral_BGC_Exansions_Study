import os
import sys

file1 = 'High_Similarity_Hits_to_Bacterial_Proteins.txt'
file2 = 'Codoff_Results_FDR_Corrected.txt'

blast_hits = set([])
with open(file1) as of1:
    for line in of1:
        line = line.strip()
        ls = line.split('\t')
        blast_hits.add(tuple([ls[0], ls[1]]))

with open(file2) as of2:
    for line in of2:
        line = line.strip()
        ls = line.split('\t')
        pair = tuple([ls[0], ls[1].split('.gbk')[0]])
        if pair in blast_hits:
            print(line)
