import os
import sys
from Bio import SeqIO
from collections import defaultdict
import multiprocessing 

contaminated = set(['JADGJL010000003.1.region001.gbk'])

def parseAntiSMASHGBKForFunctionAndCompleteness(inputs):
    gca, name, bgc, bgc_gbk, outf = inputs
    product = 'NA'
    bgc_length = 'NA'
    complete_status = True
    con_or_like = True
    try:
        with open(bgc_gbk) as obg:
            for line in obg:
                line = line.strip()
                if '/contig_edge="True"' in line:
                    complete_status = False
        products = set([])
        with open(bgc_gbk) as obg:
            for rec in SeqIO.parse(obg, 'genbank'):
                for feat in rec.features:
                    if feat.type == 'protocluster':
                        try:
                            prod = feat.qualifiers.get('product')[0]
                            if not '-like' in prod and not '-containing' in prod:
                                con_or_like = False
                            products.add(feat.qualifiers.get('product')[0])
                        except:
                            pass
                bgc_length = len(rec.seq)

        product_string = ' '.join(products)
        metallophore_related = False
        nrps_related = False
        pks_related = False
        terpene_related = False
        if 'terpene' in products:
            terpene_related = True
        if 'NRP-metallophore' in products or 'aminopolycarboxylic-acid' in products or 'opine-like-metallophore' in products or 'NI-siderophore' in products:
            metallophore_related = True
        if 'thioamide-NRP' in products or 'NRPS' in products or 'NRPS-like' in products or 'NRP-metallophore' in products:
            nrps_related = True
        if 'hgLE-KS' in products or 'PKS-like' in products or 'prodigiosin' in products or 'T1PKS' in products or 'T2PKS' in products or 'T3PKS' in products or 'transAT-PKS' in products or 'transAT-PKS-like' in products:
            pks_related = True
    except:
        sys.stderr.write('Issues parsing BGC Genbank %s\n' % bgc_gbk)
        raise RuntimeError()

    outh = open(outf, 'w')
    outh.write('\t'.join([str(x) for x in [gca, name, bgc, bgc_length, complete_status, con_or_like, metallophore_related, nrps_related, pks_related, terpene_related,  '|'.join(products), bgc_gbk]]) + '\n')
    outh.close()

asresdir = 'All_antiSMASH_Results/'
overview_file = '../Overview_File.txt'

name_mapping = {}
with open(overview_file) as ovf:
    for i, line in enumerate(ovf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        gca = ls[0]
        name = ls[1]
        name_mapping[gca] = name

listing_dir = os.path.abspath('BGC_Product_Listings/') + '/'
bgc_inputs = []
header = ['gca', 'name', 'bgc', 'bgc_length', 'completeness', 'putative', 'metallophore-like', 'nrps-like', 'pks-like', 'terpene-like', 'product-types', 'bgc_gbk']
index = 1
for i, s in enumerate(os.listdir(asresdir)):
    #if i > 3: continue
    gca = '_'.join(s.split('_')[:2])
    samp_dir = asresdir + s + '/'
    for f in os.listdir(samp_dir):
        if not f.endswith('.gbk'): continue
        if not '.region' in f: continue
        bgc_path = samp_dir + f
        if f in contaminated: continue
        name = name_mapping[gca]
        outf = listing_dir + str(index) + '.txt'
        index += 1
        bgc_inputs.append([gca, name, f.split('.gbk')[0], bgc_path, outf])

p= multiprocessing.Pool(20)
p.map(parseAntiSMASHGBKForFunctionAndCompleteness, bgc_inputs)
p.close()

print('\t'.join(header))
for f in os.listdir(listing_dir):
    with open(listing_dir + f) as olf:
        for line in olf:
            line = line.strip()
            print(line)
