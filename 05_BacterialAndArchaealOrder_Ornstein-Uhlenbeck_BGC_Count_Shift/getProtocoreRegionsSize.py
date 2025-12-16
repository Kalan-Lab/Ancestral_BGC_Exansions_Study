import os
import sys
from Bio import SeqIO
from collections import defaultdict

gcas_file = 'gcas_listing.txt'
asresdir = '/workspace/lab/kalanlab/salamzar/multicellularity_and_bgcs/researchdrive/Multicellularity_and_BGCs/antiSMASH_Order_Full_Results/'

gcas_to_investigate = set([])
with open(gcas_file) as ogf:
    for line in ogf:
        line = line.strip()
        gcas_to_investigate.add(line.split('.')[0])

gca_protocore_sums = defaultdict(lambda: 0.0)
for s in os.listdir(asresdir):
    samp_dir = asresdir + s + '/'
    gca = s.split('.')[0]
    if not gca in gcas_to_investigate: continue
    for f in os.listdir(samp_dir):
        if not f.endswith('.gbk') or '.region' not in f: continue

        core_pos = set([])
        with open(samp_dir + f) as osf:
            for rec in SeqIO.parse(osf, 'genbank'):
                for feat in rec.features:
                    if feat.type != 'proto_core': continue
                    loc = feat.location
                    start = min([loc.start, loc.end])
                    end = max([loc.start, loc.end])

                    for pos in range(start, end+1):
                        core_pos.add(pos)

        gca_protocore_sums[gca] += len(core_pos)
        
for gca in gca_protocore_sums:
    print(gca + '\t' + str(gca_protocore_sums[gca]))
