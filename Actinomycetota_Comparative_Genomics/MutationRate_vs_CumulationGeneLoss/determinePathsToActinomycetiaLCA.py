import os
import sys
from ete3 import Tree

def recursivePath(t, c, path):
    p = (t&c).up.name
    updated_path = path 
    if p != 'Node_127':
        updated_path += [p]
        updated_path = recursivePath(t, p, updated_path)
    else:
        print('\t'.join(updated_path) + '\tNode_127')
    
t = Tree(sys.argv[1], format=1)
actinomycetia = set([x.strip() for x in open(sys.argv[2]).readlines()])

for c in actinomycetia:
    path = recursivePath(t, c, [c])
