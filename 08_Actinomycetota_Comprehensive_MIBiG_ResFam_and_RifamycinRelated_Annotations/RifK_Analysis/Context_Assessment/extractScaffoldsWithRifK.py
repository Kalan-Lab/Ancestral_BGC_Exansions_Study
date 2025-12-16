import os
import sys
from Bio import SeqIO
from collections import defaultdict

rif_clade_file = 'RifK_Clade.txt'
gbk_dir = 'RifK_Search/Final_Results/Homologous_Gene_Cluster_GenBanks/'
full_gbk_dir = '/workspace/lab/kalanlab/salamzar/multicellularity_and_bgcs/researchdrive/Multicellularity_and_BGCs/All_Actino_Genomes_prepTG_DB/Genomic_Genbanks_Additional/'

rifk_genes = set([])
with open(rif_clade_file) as orcf:
    for line in orcf:
        line = line.strip()
        gca = line.split('|')[0]
        lt = line.split('|')[1]
        rifk_genes.add(lt)

gca_scaffs = defaultdict(set)
for f in os.listdir(gbk_dir):
    scaffs = set([])
    with open(gbk_dir + f) as ogf:
        for rec in SeqIO.parse(ogf, 'genbank'):
            for feat in rec.features:
                if feat.type == 'CDS':
                    lt = feat.qualifiers["locus_tag"][0] 
                    if lt in rifk_genes:
                        scaffs.add(rec.id)
                        break
    genome_gbk = full_gbk_dir + f.split('_fai')[0] + '.gbk'
    with open(genome_gbk) as oggf:
        for rec in SeqIO.parse(oggf, 'genbank'):
            if rec.id in scaffs:
                print('>' + f.split('.')[0] + '|' + rec.id + '\n' + str(rec.seq))
                break
