import os
import sys
from ete3 import Tree

gf_map_file = 'Genus_to_Family_Mapping.txt'
tree_file = 'Actino_OF_Phylogeny.tre'
family_level_tree = 'Final_Actino_HighConfidence_Phylogeny.Family_Level.tre'

family_reps = {}
rep_to_family = {}
with open(gf_map_file) as ogmf:
    for line in ogmf:
        line = line.strip()
        ls = line.split('\t')
        if ls[1] in family_reps: continue
        family_reps[ls[1]] = ls[0]
        rep_to_family[ls[0]] = ls[1]

list_to_keep = []
with open(gf_map_file) as ogmf:
    for line in ogmf:
        line = line.strip()
        ls = line.split('\t')
        if ls[0] == family_reps[ls[1]]:
            list_to_keep.append(ls[0])

t = Tree(tree_file)
t.prune(list_to_keep)

for n in t.traverse('postorder'):
    if n.is_leaf and n.name != '':
       n.name = rep_to_family[n.name] 

t.write(outfile=family_level_tree)
