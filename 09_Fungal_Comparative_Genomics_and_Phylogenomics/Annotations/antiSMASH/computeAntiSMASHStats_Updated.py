import os
import sys
from Bio import SeqIO
from collections import defaultdict

contaminated = set(['JADGJL010000003.1.region001.gbk', 'KZ989109.1.region001.gbk'])

def parseAntiSMASHGBKForFunctionAndCompleteness(bgc_gbk):
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

    return([con_or_like, metallophore_related, nrps_related, pks_related, terpene_related, bgc_length, complete_status])

asresdir = 'AntiSMASH_Results/'
overview_file = 'Overview_File.txt'

name_mapping = {}
with open(overview_file) as ovf:
    for i, line in enumerate(ovf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        gca = ls[0]
        #name = ls[1].split(';')[-1].split()[0] + '_' + gca
        name = '_'.join(ls[1].split(';')[-1].split()) + '_' + gca
        name_mapping[gca] = name

gca_complete_bgc_count = defaultdict(int)
gca_bgc_count = defaultdict(int)
gca_bgc_sum = defaultdict(int)
gca_nrp_pks_sum = defaultdict(int)
gca_ter_count = defaultdict(int)
gca_nrp_count = defaultdict(int)
gca_pks_count = defaultdict(int)
gca_met_count = defaultdict(int)
gca_total_sum = defaultdict(int)
gca_nrp_pks_count = defaultdict(int)
gca_complete_nrp_pks_count = defaultdict(int)
gca_nrp_pks_strict_sum = defaultdict(int)
for i, s in enumerate(os.listdir(asresdir)):
    #if i > 3: continue
    gca = '_'.join(s.split('_')[:2])
    samp_dir = asresdir + s + '/'
    for f in os.listdir(samp_dir):
        if not f.endswith('.gbk'): continue
        isbgc = True
        if not '.region' in f: 
            with open(samp_dir + f) as osf:
                for rec in SeqIO.parse(osf, 'genbank'):
                    gca_total_sum[gca] += len([x for x in str(rec.seq) if x.upper() in set(['A', 'C', 'G', 'T'])])
        else:
            bgc_path = samp_dir + f
            if f in contaminated: continue
            like, mr, nr, pr, tr, bl, cs = parseAntiSMASHGBKForFunctionAndCompleteness(bgc_path)
            #if like: continue
            gca_bgc_count[gca] += 1
            if cs:
                gca_complete_bgc_count[gca] += 1
            gca_bgc_sum[gca] += bl
            if mr:
                gca_met_count[gca] += 1
            if nr:
                gca_nrp_count[gca] += 1
            if pr:
                gca_pks_count[gca] += 1
            if tr:
                gca_ter_count[gca] += 1
            if nr or pr:
                gca_nrp_pks_sum[gca] += bl
                if not like:
                    gca_nrp_pks_strict_sum[gca] += bl
                gca_nrp_pks_count[gca] += 1
                if cs:
                    gca_complete_nrp_pks_count[gca] += 1

print('Assembly_ID\tAssembly_Name\tBGC_Count\tComplete_BGC_Count\tBGC_Sum\tGenome_Size\tBGC_Genome_Prop\tTerpene_Prop\tMetallophore_Prop\tNRPS_Prop\tPKS_Prop\tNRPS_or_PKS\tComplete_NRPS_or_PKS\tNRPS_or_PKS_Sum\tNRPS_or_PKS_Sum_Strict')
for gca in gca_total_sum:
    bgc_count = str(gca_bgc_count[gca])
    complete_bgc_count = str(gca_complete_bgc_count[gca])
    bgc_genome_prop = float(gca_bgc_sum[gca])/float(gca_total_sum[gca])
    
    met_prop = 0.0; nrp_prop = 0.0; pks_prop = 0.0; ter_prop = 0.0; nrp_or_pks = 0; complete_nrp_or_pks = 0
    if gca_bgc_count[gca] > 0:
        met_prop = float(gca_met_count[gca])/float(gca_bgc_count[gca])
        nrp_prop = float(gca_nrp_count[gca])/float(gca_bgc_count[gca])
        pks_prop = float(gca_pks_count[gca])/float(gca_bgc_count[gca])
        ter_prop = float(gca_ter_count[gca])/float(gca_bgc_count[gca])
        nrp_or_pks = gca_nrp_pks_count[gca]
        complete_nrp_or_pks = gca_complete_nrp_pks_count[gca]
    print('\t'.join([str(x) for x in [gca, name_mapping[gca], bgc_count, complete_bgc_count, gca_bgc_sum[gca], gca_total_sum[gca], bgc_genome_prop, ter_prop, met_prop, nrp_prop, pks_prop, nrp_or_pks, complete_nrp_or_pks, gca_nrp_pks_sum[gca], gca_nrp_pks_strict_sum[gca]]]))
