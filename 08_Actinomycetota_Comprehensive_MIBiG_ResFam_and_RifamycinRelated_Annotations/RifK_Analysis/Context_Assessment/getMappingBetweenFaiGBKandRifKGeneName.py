import os
import sys
from Bio import SeqIO
from collections import defaultdict
from ete3 import Tree

rif_clade_file = 'All_RifK_Homologs.txt'
gbk_dir = 'RifK_Search/Final_Results/Homologous_Gene_Cluster_GenBanks/'
tree_file = '../RifK_Homologs_including_MIBiG_and_Bacillota.tre'

rifk_genes = set([])
with open(rif_clade_file) as orcf:
    for line in orcf:
        line = line.strip()
        gca = line.split('|')[0]
        lt = line.split('|')[1]
        rifk_genes.add(lt)

lt_to_f = {}
for f in os.listdir(gbk_dir):
    scaffs = set([])
    with open(gbk_dir + f) as ogf:
        for rec in SeqIO.parse(ogf, 'genbank'):
            for feat in rec.features:
                if feat.type == 'CDS':
                    lt = feat.qualifiers["locus_tag"][0] 
                    start = feat.location.start
                    if lt in rifk_genes and start == 50000:
                        lt_to_f[lt] = f

t = Tree(tree_file)

for node in t.traverse('postorder'):
    if node.is_leaf:
        name = node.name
        if not '|' in name: continue
        lt = name.split('|')[1]
        if lt in lt_to_f:
            print(name + '\t' + lt_to_f[lt])
