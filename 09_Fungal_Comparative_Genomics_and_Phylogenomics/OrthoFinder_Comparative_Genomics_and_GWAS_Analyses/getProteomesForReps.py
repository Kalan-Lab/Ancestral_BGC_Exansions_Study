import os
import sys
from Bio import SeqIO

as_dir = os.path.abspath(sys.argv[1]) + '/'
res_dir = os.path.abspath(sys.argv[2]) + '/'

for s in os.listdir(as_dir):
    samp_dir = as_dir + s + '/'
    samp_prot_file = res_dir + s + '.faa'
    outf = open(samp_prot_file, 'w')
    for f in os.listdir(samp_dir):
        if not '.region' in f and f.endswith('.gbk'):
            genome_gbk = samp_dir + f
            cds_id = 1
            with open(genome_gbk) as ogf:
                for rec in SeqIO.parse(ogf, 'genbank'):
                    for feat in rec.features:
                        if feat.type != 'CDS': continue
                        translation = feat.qualifiers.get('translation')[0]
                        outf.write('>' + s + '|CDS_' + str(cds_id) + '\n' + translation + '\n')
                        cds_id += 1
    outf.close()
                                                                                                                       
