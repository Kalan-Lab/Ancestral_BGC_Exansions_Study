import os
import sys
from Bio import SeqIO
import multiprocessing 
from collections import defaultdict

genomes_dir = os.path.abspath('Genomic_Genbanks_Additional//') + '/'
gc_dir = os.path.abspath('Genomic_Temp_Estimations/') + '/'

inputs = []
for f in os.listdir(genomes_dir):
    outf = gc_dir + f + '.txt'
    inf = genomes_dir + f 
    inputs.append([inf, outf])

def mp_worker(inputs):
    inf, outf = inputs
    interest_aa = set(['I', 'V', 'Y', 'W', 'R', 'E', 'L'])

    aa_tot = 0
    aa_int = 0
    with open(inf) as osf:
        for rec in SeqIO.parse(osf, 'genbank'):
            for feature in rec.features:
                if feature.type == 'CDS':
                    aa_seq = str(feature.qualifiers.get('translation')[0]).strip('*')
                    for aa in aa_seq:
                        if aa in interest_aa: aa_int += 1
                        aa_tot += 1

    aa_frac = float(aa_int)/float(aa_tot)
    mpt = (937.0*aa_frac) - 335.0

    outf = open(outf, 'w')
    outf.write(inf + '\t' + str(aa_frac) + '\t' + str(mpt) + '\n')
    outf.close()

p = multiprocessing.Pool(40)
p.map(mp_worker, inputs)

