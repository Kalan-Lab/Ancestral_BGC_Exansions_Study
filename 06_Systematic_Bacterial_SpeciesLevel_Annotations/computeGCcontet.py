import os
import sys
from Bio import SeqIO
import multiprocessing 

genomes_dir = os.path.abspath('Genome_FASTAs/') + '/'
gc_dir = os.path.abspath('GC_Content/') + '/'

inputs = []
for f in os.listdir(genomes_dir):
    outf = gc_dir + f + '.txt'
    inf = genomes_dir + f 
    inputs.append([inf, outf])

def mp_worker(inputs):
    inf, outf = inputs
    gc = 0
    at = 0
    with open(inf) as ogf:
        for rec in SeqIO.parse(ogf, 'fasta'):
            seq = str(rec.seq).upper()
            for allele in seq:
                if allele in set(['C', 'G']):
                    gc += 1
                elif allele in set(['A', 'T']):
                    at += 1
    outf = open(outf, 'w')
    outf.write(inf + '\t' + str(float(gc)/float(gc+at)) + '\n')
    outf.close()

p = multiprocessing.Pool(40)
p.map(mp_worker, inputs)

