import os 
import sys
from Bio import SeqIO

# script to get protein sequences across genomes from antiSMASH results directory

as_resdir = os.path.abspath(sys.argv[1]) + '/'

for s in os.listdir(as_resdir):
    samp_dir = as_resdir + s + '/'
    for f in os.listdir(samp_dir):
        if f.endswith('.gbk') and not 'region' in f:
            cds_id = 1
            with open(samp_dir + f) as osf:
                for rec in SeqIO.parse(osf, 'genbank'):
                    for feature in rec.features:
                        if feature.type != 'CDS': continue
                        prot_seq = feature.qualifiers.get('translation')[0]
                        print('>' + s + '|CDS_' + str(cds_id) + '\n' + prot_seq)
                        cds_id += 1
