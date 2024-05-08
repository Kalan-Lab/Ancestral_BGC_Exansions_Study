import os
import sys
from ete3 import Tree

print('156          1')

calib_data = 'Calibration_Data.txt'

nodebounds = {}
with open(calib_data) as ocd:
    for i, line in enumerate(ocd):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        name = ls[1].strip('"')
        bound = "'B(" + ', '.join(ls[3:]) + ")'"
        nodebounds[name] = bound

print(nodebounds)
tre_file = 'Fungi_Wide_Tree.Rooted_by_Capsaspora.With_InnerNode_Names.tre'
t = Tree(tre_file, format=1)
for n in t.traverse('postorder'):
    if not n.is_leaf():
        if n.name in nodebounds:
            n.name = nodebounds[n.name]
        else:
            n.name = ''

t.write(outfile="Phylogeny_for_MCMCTree.nwk", format=1)
