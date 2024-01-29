import os
import sys
from Bio import SeqIO
from collections import defaultdict

def parseAntiSMASHGBKForFunctionAndCompleteness(bgc_gbk):
    product = 'NA'
    bgc_length = 'NA'
    complete_status = 'Not Near Contig Edge'
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

    return([metallophore_related, nrps_related, pks_related, bgc_length])

genus_file = 'Species_Level_Selections.txt'
asresdir = 'antiSMASH_Results_Phyla_Level/'

gca_to_genus = {}
with open(genus_file) as ogf:
    for line in ogf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[-1]
        gen = ls[0] + '\t' + ls[1]
        gca_to_genus[gca] = gen

gca_bgc_count = defaultdict(int)
gca_bgc_sum = defaultdict(int)
gca_nrp_sum = defaultdict(int)
gca_pks_sum = defaultdict(int)
gca_nrp_or_pks_count = defaultdict(int)
gca_met_sum = defaultdict(int)
gca_total_sum = defaultdict(int)
for s in os.listdir(asresdir):
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
            mr, nr, pr, bl = parseAntiSMASHGBKForFunctionAndCompleteness(bgc_path)
            gca_bgc_count[gca] += 1
            gca_bgc_sum[gca] += bl
            if mr:
                gca_met_sum[gca] += bl
            if nr:
                gca_nrp_sum[gca] += bl
            if pr:
                gca_pks_sum[gca] += bl
            if nr or pr:
                gca_nrp_or_pks_count[gca] += 1

print('Assembly_ID\tPhylum\tGenus\tBGC_Count\tGenome_Size\tBGC_Genome_Prop\tMetallophore_Genome_Prop\tNRPS_Genome_Prop\tPKS_Genome_Prop\tNRPS_or_PKS_Prop')
for gca in gca_total_sum:
    bgc_count = str(gca_bgc_count[gca])
    bgc_genome_prop = float(gca_bgc_sum[gca])/float(gca_total_sum[gca])
    met_genome_prop = float(gca_met_sum[gca])/float(gca_total_sum[gca])
    nrp_genome_prop = float(gca_nrp_sum[gca])/float(gca_total_sum[gca])
    pks_genome_prop = float(gca_pks_sum[gca])/float(gca_total_sum[gca])
    nrp_or_pks_genome_prop = 'NA'
    if int(bgc_count) > 0:
        nrp_or_pks_genome_prop = float(gca_nrp_or_pks_count[gca])/float(bgc_count)
    print('\t'.join([str(x) for x in [gca, gca_to_genus[gca], bgc_count, gca_total_sum[gca], bgc_genome_prop, met_genome_prop, nrp_genome_prop, pks_genome_prop, nrp_or_pks_genome_prop]]))
