import os
import sys
from ete3 import Tree

phylum = sys.argv[1]
out_file = sys.argv[2]

t = Tree('../GToTree_Results/GToTree_Results.tre') 

leaves = set([])
for n in t.traverse('postorder'):
    if n.is_leaf():
        n.name = n.name.split('.')[0]
        leaves.add(n.name)

samples_to_keep = set([])
with open('Comprehensive_Info_for_Bacterial_Species_Reps.tsv') as obs:
    for line in obs:
        line = line.strip()
        ls = line.split('\t')
        if ls[0] in leaves and ls[5] == phylum:
           samples_to_keep.add(ls[0])

print(len(samples_to_keep))

t.prune(samples_to_keep, preserve_branch_length=True)
t.write(outfile=out_file)
