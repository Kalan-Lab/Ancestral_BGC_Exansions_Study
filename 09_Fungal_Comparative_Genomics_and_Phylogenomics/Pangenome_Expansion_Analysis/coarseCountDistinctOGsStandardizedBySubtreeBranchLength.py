import os
import sys
from collections import defaultdict
from Bio import SeqIO
from ete3 import Tree
import statistics
from operator import itemgetter
import random

og_tsv_file = 'Orthogroups.tsv'
sog_tsv_file = 'Orthogroups_UnassignedGenes.tsv'
clades_dir = 'Clades/'
subtree_dir = 'Subtrees/'

gca_ogs = defaultdict(set)
genomes = []
with open(og_tsv_file) as ootf:
    for i, line in enumerate(ootf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            genomes = ls[1:]
        else:
            og = ls[0]
            for j, gs in enumerate(ls[1:]):
                g = genomes[j]
                if gs.strip() != '':
                    gca_ogs[g].add(og)

genomes = []
with open(sog_tsv_file) as ootf:
    for i, line in enumerate(ootf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            genomes = ls[1:]
        else:
            og = ls[0]
            for j, gs in enumerate(ls[1:]):
                g = genomes[j]
                if gs.strip() != '':
                    gca_ogs[g].add(og)

taxa_gcas = defaultdict(set)
for f in os.listdir(clades_dir):
    taxa_file = clades_dir + f
    taxa = f.split('.txt')[0] 
    with open(taxa_file) as otf:
        for line in otf:
            line = line.strip()
            if '_GCA' in line:
                gca = 'GCA_' + line.split('_GCA_')[1]
                #if gca in reps:
                taxa_gcas[taxa].add(gca)

select_taxa = ['All_Fungi_Genomes', 'NonYeast_Basidi', 'Zoopago', 'Mucoro', 'Spizello', 'Neocalli', 'Pezizomycotina', 'Zoosporic_or_Basal', 'Yeast_Basidi', 'Yeast_Ascomyco', 'Chytridiales', 'Rhizophydiales', 'Blasto', 'Zygomycota', 'Dikarya']
select_taxa_name = ['All Fungal Genomes', 'Non-yeast Baidiomycota', 'Zoopagomycota', 'Mucoromycota', 'Spizellomycetales & Rhizophlyctidales', 'Neocallimastigomycota', 'Pezizomycotina', 'Zoosporic or Basal Fungi', 'Yeast Basidiomycota', 'Yeast Ascomycota', 'Chytridiales', 'Rhizophydiales', 'Blastocladiomycota', 'Zygomycota', 'Dikarya']
select_taxa_type = ['Variable', 'Diploid+ Dominant', 'Diploid+ Dominant', 'Haploid Dominant', 'Haploid Dominant', 'Haploid Dominant', 'Haploid Dominant', 'Variable', 'Variable', 'Variable', 'Diploid+ Dominant', 'Diploid+ Dominant', 'Diploid+ Dominant', 'Variable', 'Variable']

print('\t'.join(['statistic', 'clade', 'clade_ploidy', 'number_of_genomes', 'avg_ogs', 'total_ogs', 'total_phylo_distance']))
for j, st in enumerate(select_taxa):
    gcas = list(taxa_gcas[st])
    total_ogs = set([])
    og_counts = []
    for gca in sorted(gcas):
        total_ogs = total_ogs.union(gca_ogs[gca])
        og_counts.append(len(gca_ogs[gca]))

    phylo_breadth = 0.0
    st = Tree(subtree_dir + st + '.tre')

    for n in st.traverse('postorder'):
        phylo_breadth += n.dist

    avg_og_count = statistics.mean(og_counts)
    og_stand = (len(total_ogs)-avg_og_count)/phylo_breadth
    print('\t'.join([str(x) for x in [round(og_stand, 3), select_taxa_name[j], select_taxa_type[j], len(gcas), avg_og_count, len(total_ogs), phylo_breadth]]))
