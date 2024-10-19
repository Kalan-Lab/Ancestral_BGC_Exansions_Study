import os
import sys
from ete3 import Tree

gca_to_name = {}
with open('Overview_File.txt') as ovf:
    for line in ovf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        gca_to_name[gca] = ls[1]

t = Tree('IQTree_SCC.treefile')
for n in t.traverse('postorder'):
    if n.is_leaf():
        gca = n.name
        sp_name = gca_to_name[gca]
        n.name = sp_name

t.write(outfile='Fungi_Wide_Tree_with_Outgroups.tre')
