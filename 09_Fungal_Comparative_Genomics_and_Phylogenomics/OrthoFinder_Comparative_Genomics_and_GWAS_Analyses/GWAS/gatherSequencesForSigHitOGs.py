import os
import sys
from Bio import SeqIO

sig_hits = set([])
with open('Significant_Hits.txt') as oshf:
    for line in oshf:
        line = line.strip()
        ls = line.split('\t')
        sig_hits.add(ls[0])

og_seq_dir = '../../OrthoFinder_Analysis/OrthoFinder_Results/Results_Oct02/Orthogroup_Sequences/'

for f in os.listdir(og_seq_dir):
    og = f.split('.fa')[0]
    if not og in sig_hits: continue
    with open(og_seq_dir + f) as osf:
        for rec in SeqIO.parse(osf, 'fasta'):
            print('>' + og + '|' + rec.id + '\n' + str(rec.seq))
