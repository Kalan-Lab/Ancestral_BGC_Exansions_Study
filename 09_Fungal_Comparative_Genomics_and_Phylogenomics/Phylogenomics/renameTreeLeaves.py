import os
import sys
from ete3 import Tree

naming = {}
with open('Overview_File.txt') as ovf:
    for line in ovf:
        line = line.strip()
        ls = line.split('\t')
        naming[ls[0]] = ls[1]

t = Tree("IQTRee_GHOST_SCC.treefile")
for n in t.traverse('postorder'):
    if n.is_leaf():
        n.name = naming[n.name]

t.write(outfile="IQTree_Renamed.tre", format=1)

