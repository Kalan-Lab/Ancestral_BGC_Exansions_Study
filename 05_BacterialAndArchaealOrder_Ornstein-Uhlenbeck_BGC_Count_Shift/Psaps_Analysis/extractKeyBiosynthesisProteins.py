import os
import sys
from Bio import SeqIO
from collections import defaultdict

SEARCH_STRING = "biosynthetic (rule-based-clusters"

gcas_file = 'select_gcas_without_issues.txt'
asresdir = '/workspace/lab/kalanlab/salamzar/multicellularity_and_bgcs/researchdrive/Multicellularity_and_BGCs/antiSMASH_Order_Full_Results/'

gcas_to_investigate = set([])
with open(gcas_file) as ogf:
    for line in ogf:
        line = line.strip()
        gcas_to_investigate.add(line.split('.')[0])

for s in os.listdir(asresdir):
    samp_dir = asresdir + s + '/'
    gca = s.split('.')[0]
    if not gca in gcas_to_investigate: continue
    for f in os.listdir(samp_dir):
        if not f.endswith('.gbk') or '.region' not in f: continue

        with open(samp_dir + f) as osf:
            for rec in SeqIO.parse(osf, 'genbank'):
                for feat in rec.features:
                    if feat.type != 'CDS': continue
                    lt = feat.qualifiers.get("locus_tag")[0]

                    if "gene_functions" in feat.qualifiers:
                        function_list = feat.qualifiers["gene_functions"]

                        has_biosynthetic_func = any(func.replace("\n", " ").strip().startswith(SEARCH_STRING) 
                                                for func in function_list)

                        if has_biosynthetic_func:
                            prot_seq = feat.qualifiers.get("translation")[0]
                            print('>' + gca + '|' + lt + '\n' + prot_seq)
                    
