import os
import sys
from ete3 import Tree

t = Tree('IQTree_SCC.treefile')
selection = [x.strip() for x in open('Select_Genomes_for_TimeTree_55markers_or_more.txt').readlines()]
print(selection)
t.prune(selection, preserve_branch_length=True)
t.set_outgroup('GCA_000151315')

naming = {}
with open('Overview_File.txt') as ovf:
    for line in ovf:
        line = line.strip()
        ls = line.split('\t')
        naming[ls[0]] = ls[1]

for n in t.traverse('postorder'):
    if n.is_leaf():
        gca = n.name
        n.name= naming[gca]

t.write(outfile='TimeTree_Selection_PBL.tre', format=1)
