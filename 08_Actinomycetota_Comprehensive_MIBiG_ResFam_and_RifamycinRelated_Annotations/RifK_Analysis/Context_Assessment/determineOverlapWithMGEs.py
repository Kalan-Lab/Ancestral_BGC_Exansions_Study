import os
import sys
from collections import defaultdict
from Bio import SeqIO

full_gbk_dir = '/workspace/lab/kalanlab/salamzar/multicellularity_and_bgcs/researchdrive/Multicellularity_and_BGCs/All_Actino_Genomes_prepTG_DB/Genomic_Genbanks_Additional/'
rifk_file = 'RifK_to_FaiGC.txt'
virus_file = 'geNomad_Results/Scaffolds_with_RifK_Homologs_summary/Scaffolds_with_RifK_Homologs_virus_summary.tsv'
plasmid_file = 'geNomad_Results/Scaffolds_with_RifK_Homologs_summary/Scaffolds_with_RifK_Homologs_plasmid_summary.tsv'

gca_to_gbk = {}
for f in os.listdir(full_gbk_dir): 
    gca = f.split('.')[0]
    gca_to_gbk[gca] = full_gbk_dir + f

gca_mge_scaffs = defaultdict(set)
gca_mge_scaff_coords = defaultdict(lambda: defaultdict(set))
with open(virus_file) as ovf:
    for i, line in enumerate(ovf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t') 
        if 'provirus_' in ls[0]:
            gca, scaff, coords = ls[0].split('|')
            start = int(coords.split('_')[-2])
            end = int(coords.split('_')[-1])
            pp_range = tuple([start, end])
            gca_mge_scaff_coords[gca][scaff].add(pp_range)
        else:
            gca, scaff = ls[0].split('|')
            gca_mge_scaffs[gca].add(scaff)

gca_plas_scaffs = defaultdict(set)
with open(plasmid_file) as opf:
    for i, line in enumerate(opf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca, scaff = ls[0].split('|')
        gca_plas_scaffs[gca].add(scaff)

print('DATASET_COLORSTRIP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tgeNomad')
print('DATA')
with open(rifk_file) as orf:
    for line in orf:
        line = line.strip()
        ls = line.split('\t')
        gc = ls[1].split('.gbk')[0]
        gca = gc.split('.')[0]
        gbk = gca_to_gbk[gca]
        lt = ls[0].split('|')[1]
        with open(gbk) as ogbk:
            for rec in SeqIO.parse(ogbk, 'genbank'):
                scaff = rec.id
                if scaff in gca_mge_scaffs[gca] or scaff in gca_mge_scaff_coords[gca]:
                    for feat in rec.features:
                        if feat.type == 'CDS':
                            cds_lt = feat.qualifiers['locus_tag'][0]
                            if cds_lt == lt:
                                if scaff in gca_plas_scaffs[gca]:
                                    print(ls[0] + '\t#c25f64\tplasmid')
                                elif scaff in gca_mge_scaffs[gca]:
                                    #print(gca)
                                    #print(scaff)
                                    #print(feat.location)
                                    print(ls[0] + '\t#c25f6450\tphage')
                                    #print('---------------')
                                else:
                                    cds_coords = set([])
                                    for position in feat.location:
                                        cds_coords.add(position+1)
                                    for pp in gca_mge_scaff_coords[gca][scaff]:
                                        pp = list(pp)
                                        pp_coords = set(range(pp[0], pp[1]+1))
                                        if len(cds_coords.intersection(pp_coords)) > 0:
                                            #print(gca)
                                            #print(scaff)
                                            #print(feat.location)
                                            #print(pp)
                                            print(ls[0] + '\t#c25f6450\tprophage')
                                            #print('-------')
