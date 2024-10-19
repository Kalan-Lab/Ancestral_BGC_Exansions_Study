import os
import sys
from ete3 import Tree

t = Tree(sys.argv[1])
g1 = sys.argv[2]
g2 = sys.argv[3]

print(t.get_distance(g1, g2))

t.prune([g1, g2], preserve_branch_length=True)
phylo_breadth = 0.0
print('-------------------')
for n in t.traverse():
    if n.is_root(): continue
    phylo_breadth += n.dist
print('----------------')
print(phylo_breadth)


print(t.get_distance(g1, g2))
