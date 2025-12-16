import os
import sys
from Bio import SeqIO

gbk_dir = 'RifK_Search/Final_Results/Homologous_Gene_Cluster_GenBanks/'
keep_dir = 'RifK_Context_GBKs/'
rif_list_file = 'All_RifK_Homologs.txt'

rifk_genes = set([])
with open(rif_list_file) as orcf:
    for line in orcf:
        line = line.strip()
        gca = line.split('|')[0]
        lt = line.split('|')[1]
        rifk_genes.add(lt)

found = set([])
for f in os.listdir(gbk_dir):
    match = False
    with open(gbk_dir + f) as ogf:
        for rec in SeqIO.parse(ogf, 'genbank'):
            seqlen = len(rec.seq)
            if seqlen < 100000: continue
            for feat in rec.features:
                if feat.type == 'CDS':
                    lt = feat.qualifiers["locus_tag"][0] 
                    start = feat.location.start
                    if start == 50000 and lt in rifk_genes:
                        match = True
                        found.add(lt)
                        break
    if match:
        os.system('cp ' + gbk_dir + f + ' ' + keep_dir)

for g in rifk_genes:
    if not g in found:
        print(g)
