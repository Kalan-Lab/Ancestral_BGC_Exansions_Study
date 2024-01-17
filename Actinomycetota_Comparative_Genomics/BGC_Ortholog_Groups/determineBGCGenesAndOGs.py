import os
import sys
from collections import defaultdict
from Bio import SeqIO

def parseAntiSMASHGBKForProductInfo(bgc_gbk):
    products_by_location = defaultdict(set)
    try:
        with open(bgc_gbk) as obg:
            for rec in SeqIO.parse(obg, 'genbank'):
                for feat in rec.features:
                    if feat.type == 'protocluster':
                        start = min([int(x.strip('>').strip('<')) for x in str(feat.location)[1:].split(']')[0].split(':')]) + 1
                        end = max([int(x.strip('>').strip('<')) for x in str(feat.location)[1:].split(']')[0].split(':')])
                        try:
                            product = feat.qualifiers.get('product')[0]
                            for pos in range(start, end+1):
                                products_by_location[pos].add(product)
                        except:
                            pass
    except:
        print('error processing genbank')
        sys.exit(1)
    return(products_by_location)

og_mat_file = '../OrthoFinder_Results/Orthogroups.tsv'
common_file = '../Representative_Genomes_for_Comparative_Genomics.txt'
proteomes_og_dir = '../Proteomes/' # Not included - but it is a directory featuring CDS sequences for representative genomes in FASTA format 
as_dir = 'antiSMASH_Results/' # note, results compressed in tar.gz to reduce size

lt_to_og = {}
with open(og_mat_file) as omf:
    for i, line in enumerate(omf):
        if i == 0: continue
        line = line.strip('\n')
        ls = line.split('\t')
        for gs in ls[1:]:
            for g in gs.split(', '):
                glt = g.split('|')[-1]
                lt_to_og[glt] = ls[0]

prot_to_seq = defaultdict(dict)
prot_gca_to_lt = {}
for f in os.listdir(proteomes_og_dir):
    prot_file = proteomes_og_dir + f 
    gca = f.split('.')[0]
    with open(prot_file) as opf:
        for rec in SeqIO.parse(opf, 'fasta'):
            prot_to_seq[gca][rec.id.split('|')[-1]] = str(rec.seq)
            prot_gca_to_lt[gca] = rec.id.split('|')[-1].split('_')[0]

gca_to_genus = {}
with open(common_file) as ocf:
    for line in ocf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[-1].split('.')[0]
        genus = ls[0]
        gca_to_genus[gca] = genus

for s in os.listdir(as_dir):
    gca = s.split('.')[0]
    if gca in gca_to_genus and gca in prot_gca_to_lt:
        genus = gca_to_genus[gca]
        samp_dir = as_dir + s + '/'
        for f in os.listdir(samp_dir):
            if not f.endswith('.gbk') or not 'region' in f: continue
            gbk_file = samp_dir + f            

            products_by_location = parseAntiSMASHGBKForProductInfo(gbk_file)

            with open(gbk_file) as ogf:
                for rec in SeqIO.parse(ogf, 'genbank'):
                    for feat in rec.features:
                        if feat.type != 'CDS': continue
                        lt = feat.qualifiers.get('locus_tag')[0]
                        rule_based_bgc_cds = False
                        try:
                            if 'rule-based-clusters' in feat.qualifiers.get('gene_functions')[0]:
                                rule_based_bgc_cds = True
                        except:
                            pass
                        if not rule_based_bgc_cds: continue
                        olt = prot_gca_to_lt[gca] + '_' + lt.split('_')[1]
                        translation = feat.qualifiers.get('translation')[0]
                        start = min([int(x.strip('>').strip('<')) for x in str(feat.location)[1:].split(']')[0].split(':')]) + 1
                        end = max([int(x.strip('>').strip('<')) for x in str(feat.location)[1:].split(']')[0].split(':')])
                        products = set([])
                        for pos in range(start, end+1):
                            products = products.union(products_by_location[pos])
                        product = ' '.join(products)
                        if olt in lt_to_og:
                            assert(translation == prot_to_seq[gca][olt])
                            print(gca + '\t' + genus + '\t' + olt + '\t' + lt_to_og[olt] + '\t' + f.split('.gbk')[0] + '\t' + product)
                        else:
                            print(gca + '\t' + genus + '\t' + olt + '\tNA\t' + f.split('.gbk')[0] + '\t' + product)
