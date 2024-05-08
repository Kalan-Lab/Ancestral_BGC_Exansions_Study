import os
import sys
from Bio import SeqIO

def parseAntiSMASHGBKForFunctionAndCompleteness(bgc_gbk):
    product = 'NA'
    products = set([])
    with open(bgc_gbk) as obg:
        for rec in SeqIO.parse(obg, 'genbank'):
            for feature in rec.features:
                if feature.type == 'protocluster':
                    try:
                        products.add(feature.qualifiers.get('product')[0])
                    except:
                        pass

    terpene_related = False
    if 'terpene' in products:
        terpene_related = True
    return([terpene_related, '; '.join(products)])

overview_file = '../../Overview_File.txt'
as_dir = os.path.abspath('../../AntiSMASH_Results/') + '/'

gca_to_genus = {}
with open(overview_file) as ovf:
    for line in ovf:
        line = line.strip()
        ls = line.split('\t')
        gca_to_genus[ls[0]] = ls[1].split(';')[-1].split()[0]

contam = set(['JADGJL010000003.1.region001.gbk', 'KZ989109.1.region001.gbk']) 

count = 1
for s in os.listdir(as_dir):
    samp_dir = as_dir + s + '/'
    for f in os.listdir(samp_dir):
        if '.region' in f and f.endswith('.gbk'):
            gbk_file = samp_dir + f
            if f in contam: continue
            terpene_related, product_string = parseAntiSMASHGBKForFunctionAndCompleteness(gbk_file)
            if terpene_related:
                with open(gbk_file) as ogf:
                    for rec in SeqIO.parse(ogf, 'genbank'):
                        for feature in rec.features:
                            if feature.type == 'CDS':
                                sequence = feature.qualifiers["translation"][0]
                                label = 'NA'
                                try:
                                    label = feature.qualifiers["locus_tag"][0]
                                except:
                                    label = feature.qualifiers["gene"][0]

                                rule_based_bgc_cds = False
                                try:
                                    if 'rule-based-clusters' in feature.qualifiers.get('gene_functions')[0]:
                                        rule_based_bgc_cds = True
                                except:
                                    pass


                                if rule_based_bgc_cds:
                                    print('>' + str(count) + ' ' + gca_to_genus[s] + '|' + s + '|' + label + '\n' + sequence)
                                    count += 1
