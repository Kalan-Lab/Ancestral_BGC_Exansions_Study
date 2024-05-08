import os
import sys
from ete3 import Tree

count_file = 'AW_g2.txt'
tree_file = 'RubrobacterD_Rooted_Tree.Innernodes_Named.tre'
path_to_lca_file = 'Paths_of_Actinomycetia_to_LCA.txt'
taxa_file = 'GCA_to_Taxonomy.txt'

c1 = set([])
c2 = set([])

with open('Clade-1.txt') as oc1:
    for line in oc1:
        line = line.strip()
        c1.add(line)

with open('Clade-2.txt') as oc2:
    for line in oc2:
        line = line.strip()
        c2.add(line)

genus_fam = {}
with open(taxa_file) as otf:
    for line in otf:
        line = line.strip()
        taxa = line.split('\t')[1]
        genus = taxa.split(';g__')[1].split(';s__')[0]
        fam = taxa.split(';f__')[1].split(';g__')[0]
        genus_fam[genus] = fam

t = Tree(tree_file, format=1)

node_loss = {}
with open(count_file) as ocf:
    for line in ocf:
        line = line.strip()
        ls = line.split('\t')
        if ls[0] == '# CHANGE':
            if ls[1] == 'node' or ls[1] == 'total': continue
            node = ls[1]
            gene_loss = ls[6]
            node_loss[node] = int(ls[6])
            
with open(path_to_lca_file) as oplf:
    for i, line in enumerate(oplf):
        line = line.strip()
        ls = line.split('\t')
        genus = ls[0]
        sum_loss = 0
        for n in ls[:-1]:
            sum_loss += node_loss[' '.join(n.split('_'))]
        genus_to_lca = t.get_distance(genus, 'Node_127')
        clade = 'NA'
        if genus in c1:
            clade = 'Clade-1'
        elif genus in c2:
            clade = 'Clade-2'
        print(genus + '\t' + genus_fam[genus] + '\t' + clade + '\t' + str(genus_to_lca) + '\t' + str(sum_loss))
