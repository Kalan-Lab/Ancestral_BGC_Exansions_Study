import os
import sys
from Bio import SeqIO

mibig_gbk_dir = os.path.abspath(sys.argv[1]) + '/'

for f in os.listdir(mibig_gbk_dir):
    gbk_file = mibig_gbk_dir + f
    with open(gbk_file) as ogf:
        for rec in SeqIO.parse(gbk_file, 'genbank'):
            for feature in rec.features:
                if feature.type != 'CDS': continue
                gene_type = None
                try:
                    gene_type = feature.qualifiers.get('gene_kind')[0]
                except:
                    continue
                if gene_type != 'biosynthetic': continue
                lt = None
                try:
                    lt = feature.qualifiers.get('locus_tag')[0]
                except:
                    lt = feature.qualifiers.get('protein_id')[0]
                assert(lt != None)
                prot_seq = feature.qualifiers.get('translation')[0]
                product = "NA"
                try:
                    product = feature.qualifiers.get('product')[0]
                except:
                    pass
                print('>' + f.split('.gbk')[0] + '|' + lt + '|' + product + '\n' + str(prot_seq))
