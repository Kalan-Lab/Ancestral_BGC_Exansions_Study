import os
import sys
from Bio import SeqIO
from collections import defaultdict

def parseAntiSMASHGBKForFunctionAndCompleteness(bgc_gbk):
    product = 'NA'
    bgc_length = 'NA'
    complete_status = 'Not Near Contig Edge'
    like_or_con = True
    try:
        with open(bgc_gbk) as obg:
            for line in obg:
                line = line.strip()
                if '/contig_edge="True"' in line:
                    complete_status = 'Near Contig Edge'
        products = set([])
        with open(bgc_gbk) as obg:
            for rec in SeqIO.parse(obg, 'genbank'):
                for feat in rec.features:
                    if feat.type == 'protocluster':
                        try:
                            prod = feat.qualifiers.get('product')[0]
                            if not '-like' in prod and not '-containing' in prod:
                                like_or_con = False
                            products.add(feat.qualifiers.get('product')[0])
                        except:
                            pass
                bgc_length = len(rec.seq)

        product_string = ' '.join(products)
        metallophore_related = False
        nrps_related = False
        pks_related = False
        if 'NRP-metallophore' in products or 'aminopolycarboxylic-acid' in products or 'opine-like-metallophore' in products or 'NI-siderophore' in products:
            metallophore_related = True
        if 'thioamide-NRP' in products or 'NRPS' in products or 'NRPS-like' in products or 'NRP-metallophore' in products:
            nrps_related = True
        if 'hgLE-KS' in products or 'PKS-like' in products or 'prodigiosin' in products or 'T1PKS' in products or 'T2PKS' in products or 'T3PKS' in products or 'transAT-PKS' in products or 'transAT-PKS-like' in products:
            pks_related = True
    except:
        sys.stderr.write('Issues parsing BGC Genbank %s\n' % bgc_gbk)
        raise RuntimeError()

    return([like_or_con, metallophore_related, nrps_related, pks_related, bgc_length])

genus_file = 'Selected_Genomes.txt'
asresdir = '/researchdrive_files/Personal/Rauf/Multicellularity_and_BGCs/antiSMASH_Order_Full_Results/'

gca_to_taxa = {}
with open(genus_file) as ogf:
    for line in ogf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[-1] 
        order = ls[0]
        genus = ls[1]
        gca_to_taxa[gca] = order + '\t' + genus

gca_bgc_count = defaultdict(int)
gca_bgc_sum = defaultdict(int)
gca_nrp_sum = defaultdict(int)
gca_pks_sum = defaultdict(int)
gca_nrp_or_pks_sum = defaultdict(int)
gca_met_sum = defaultdict(int)
gca_total_sum = defaultdict(int)
for i, s in enumerate(os.listdir(asresdir)):
    #if i > 10: continue
    gca = '_'.join(s.split('_')[:2])
    samp_dir = asresdir + s + '/'
    for f in os.listdir(samp_dir):
        isbgc = True
        if not '.region' in f: 
            with open(samp_dir + f) as osf:
                for rec in SeqIO.parse(osf, 'genbank'):
                    gca_total_sum[gca] += len(str(rec.seq))
        else:
            bgc_path = samp_dir + f
            like, mr, nr, pr, bl = parseAntiSMASHGBKForFunctionAndCompleteness(bgc_path)
            if like: continue
            gca_bgc_count[gca] += 1
            gca_bgc_sum[gca] += bl
            if mr:
                gca_met_sum[gca] += bl
            if nr:
                gca_nrp_sum[gca] += bl
            if pr:
                gca_pks_sum[gca] += bl
            if nr or pr:
                gca_nrp_or_pks_sum[gca] += bl

print('Assembly_ID\tOrder\tGenus\tBGC_Count\tBGCome_Size\tGenome_Size\tBGC_Genome_Prop\tMetallophore_Size\tNRPS_Size\tPKS_Size\tNRPS_or_PKS_Size')
for gca in gca_total_sum:
    bgc_count = str(gca_bgc_count[gca])
    bgc_genome_prop = float(gca_bgc_sum[gca])/float(gca_total_sum[gca])
    print('\t'.join([str(x) for x in [gca, gca_to_taxa[gca], bgc_count, gca_bgc_sum[gca], gca_total_sum[gca], bgc_genome_prop, gca_met_sum[gca], gca_nrp_sum[gca], gca_pks_sum[gca], gca_nrp_or_pks_sum[gca]]]))
