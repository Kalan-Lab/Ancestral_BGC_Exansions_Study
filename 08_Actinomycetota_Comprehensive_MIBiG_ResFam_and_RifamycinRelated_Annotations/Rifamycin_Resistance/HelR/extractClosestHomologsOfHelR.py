import os
import sys
from Bio import SeqIO
from collections import defaultdict

blast_results_file = 'Search_Full_Actinos_for_Rifamycin_Resistance_Prots.txt'
gca_to_tax_file = '../GCA_to_Taxonomy.txt'
all_prots_fasta_file = 'All_Actinos_in_GTDB-R214_Proteins.faa'
species_reps_file = 'Species_Level_Selections.GCA_IDs.txt'

species_reps = set([])
with open(species_reps_file) as osrf:
    for line in osrf:
        line = line.strip()
        species_reps.add(line.split('.')[0])

gca_to_fam = {}
fam_to_cla = {}
with open(gca_to_tax_file) as ogttf:
    for line in ogttf:
        line = line.strip()
        ls = line.split('\t')
        fam = ls[1].split('f__')[1].split(';g__')[0]
        cla = ls[1].split('c__')[1].split(';o__')[0]
        gca_to_fam[ls[0]] = fam
        fam_to_cla[fam] = cla

selected_prots = set([])
with open(blast_results_file) as obrf:
    for line in obrf:
        line = line.strip()
        ls = line.split('\t')
        if ls[0] == 'gb|CCA59315.1|-|HelR':
            qcov = float(ls[-1])
            scov = float(ls[-2])
            pid = float(ls[2])
            if qcov >= 50.0 and scov >= 50.0 and pid >= 30.0:
                bs = float(ls[-3])
                hit = ls[1]
                gca = '_'.join(ls[1].split('_')[:2]).split('.')[0]
                if not gca in species_reps: continue
                selected_prots.add(hit)

with open(all_prots_fasta_file) as oapff:
    for rec in SeqIO.parse(oapff, 'fasta'):
        if rec.id in selected_prots:
            gca = '_'.join(rec.id.split('_')[:2]).split('.')[0]
            fam = gca_to_fam[gca]
            clas = fam_to_cla[fam]
            print('>' + gca_to_fam[gca] + '|' + rec.id + '\n' + str(rec.seq))
