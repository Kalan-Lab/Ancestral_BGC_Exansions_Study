import os
import sys
from Bio import SeqIO

def parseAntiSMASHGBKForFunctionAndCompleteness(bgc_gbk):
    product = 'NA'
    products = set([])
    with open(bgc_gbk) as obg:
        for rec in SeqIO.parse(obg, 'genbank'):
            for feat in rec.features:
                if feat.type == 'protocluster':
                    try:
                        products.add(feat.qualifiers.get('product')[0])
                    except:
                        pass
        
        products_filt = products.intersection(set(['hgLE-KS', 'PKS-like', 'prodigiosin', 'T1PKS', 'T2PKS', 'T3PKS', 'transAT-PKS', 'transAT-PKS-like']))
        product_string = ';'.join(sorted(products_filt))
        pks_related = False
        if 'hgLE-KS' in products or 'PKS-like' in products or 'prodigiosin' in products or 'T1PKS' in products or 'T2PKS' in products or 'T3PKS' in products or 'transAT-PKS' in products or 'transAT-PKS-like' in products:
            pks_related = True
        return([pks_related, product_string])

as_dir = 'AntiSMASH_Results/'

for s in os.listdir(as_dir):
    samp_dir = as_dir + s + '/'
    for f in os.listdir(samp_dir):
        if '.region' in f and f.endswith('.gbk'):
            gbk_file = samp_dir + f
            pks_related, product_string = parseAntiSMASHGBKForFunctionAndCompleteness(gbk_file)

            if pks_related:
                with open(gbk_file) as ogf:
                    for rec in SeqIO.parse(ogf, 'genbank'):
                        for feat in rec.features:
                            if feat.type == 'aSDomain':
                                domain_type = feat.qualifiers["aSDomain"][0]
                                if '_KS' in domain_type:#domain_type == 'PKS_KS':
                                    label = feat.qualifiers["label"][0]
                                    sequence = feat.qualifiers["translation"][0] 
                                    print('>' + s + '|' + f.split('.gbk')[0] + '|' + product_string + '|' + label + '\n' + str(sequence))
