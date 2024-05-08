import os
import sys
from collections import defaultdict

aafreqs_file = 'Fungal_Amino_Acid_Frequencies.txt'
gc_file = 'Fungal_Genome_GC_Contents.txt'
as_stats_file = 'AntiSMASH_Stats.Updated_WithoutContaminants.txt'
cazy_file = 'CAZy_Total_Homologs.iTol.txt'

pez_clade_file = 'Clades/Pezizomycotina.txt'
yea_clade_file = 'Clades/Dikarya_Yeasts.txt'
bas_clade_file = 'Clades/NonYeast_Basidi.txt'
neo_clade_file = 'Clades/Neocalli.txt'
zoo_clade_file = 'Clades/Zoopago.txt'
muc_clade_file = 'Clades/Mucoro.txt'
spi_clade_file = 'Clades/Spizello.txt'
out_clade_file = 'Clades/Other_Zoosporic_or_Basal.txt'

genus_reps_file = 'Genus_Representative_GCAs.txt'

general_info_file = 'General_Info_on_Assemblies_and_Species.txt'

genus_reps = set([])
with open(genus_reps_file) as ogrf:
    for line in ogrf:
        line = line.strip()
        genus_reps.add(line)


clade_lists = [pez_clade_file, bas_clade_file, yea_clade_file, muc_clade_file, zoo_clade_file, neo_clade_file, spi_clade_file, out_clade_file]
clade_names = ['Pezizomycotina', 'Non-Yeast Basidimycota', 'Yeast Dikarya', 'Mucoromycota', 'Zoopagomycota', 'Neocallimastigomycota', 'Spizellomycetales', 'Other Zoosporic and Basal Fungi']

clade_gcas = defaultdict(set)
gca_to_clade = defaultdict(lambda: 'NA')
for i, clf in enumerate(clade_lists):
    cn = clade_names[i]
    with open(clf) as ocf:
        for line in ocf:
            line = line.strip()
            gca = '_'.join(line.split('_')[-2:])
            clade_gcas[cn].add(gca)
            gca_to_clade[gca] = cn

general_info = {}
with open(general_info_file) as ogif:
    for line in ogif:
        line = line.strip('\n')
        ls = line.split('\t')
        general_info[ls[0].split('.')[0].replace('GCF_', 'GCA_')] = ls[1:]

cazy_counts = defaultdict(int)
with open(cazy_file) as ocf:
    for i, line in enumerate(ocf):
        line = line.strip()
        if i > 5:
            ls = line.split('\t')
            gca = '_'.join(ls[0].split('_')[-2:])
            cazy_counts[gca] = int(ls[1])

aa_freqs = {}
with open(aafreqs_file) as oaf:
    for i, line in enumerate(oaf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        aa_freqs[ls[0]] = ls[1:]

gcs = {}
with open(gc_file) as ogf:
    for line in ogf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        gcs[gca] = ls[1]

print('\t'.join(['gca', 'name', 'Selection basis', 'Species', 'Full taxonomy', 'Annotation source/method', 'clade', 'is_genus_rep', 'bgcome_size', 'strict_nrps_and_pks_ome_size', 'genome_size', 'gc', 'cazy_counts', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])) 

with open(as_stats_file) as ogf:
    for i, line in enumerate(ogf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        name = ls[1]
        genome_size = ls[5]
        bgcome_size = ls[4]
        strict_nrps_and_pks_ome_size = ls[-1]
        clade = gca_to_clade[gca]

        cazy_count = cazy_counts[gca]
        gc = gcs[gca]
        is_genus_rep = False
        if gca in genus_reps:
            is_genus_rep = True
        print('\t'.join([str(x) for x in ([gca, name] + general_info[gca] + [clade, is_genus_rep, bgcome_size, strict_nrps_and_pks_ome_size, genome_size, gc, cazy_count] + aa_freqs[gca])]))
