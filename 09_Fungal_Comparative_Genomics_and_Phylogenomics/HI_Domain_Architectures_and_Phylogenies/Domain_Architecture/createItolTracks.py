import os
import sys
from collections import defaultdict
from ete3 import Tree

dom_list_file = sys.argv[1]
select_dom_file = sys.argv[2]
tree_file = sys.argv[3] 

select_doms = []
with open(select_dom_file) as osdf:
    for line in osdf:
        if line.startswith('standalone'): continue
        line = line.strip()
        ls = line.split()
        select_doms.append(ls[0])

prot_doms = defaultdict(set)
with open(dom_list_file) as odlf:
    for line in odlf:
        line = line.strip()
        ls = line.split('\t')
        prot = ls[0]#.split('|')[0]
        for d in ls[1].split('|'):
            prot_doms[prot].add(d)

print('DATASET_HEATMAP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\t' + sys.argv[4])
print('FIELD_LABELS\t' + '\t'.join(select_doms))
print('DATA')

t = Tree(tree_file)
for n in t.traverse('postorder'):
    if n.is_leaf():
        prot = '|'.join(n.name.split('|')[:2])
        printlist = [n.name]
        for sd in select_doms:
            if sd in prot_doms[prot]: printlist.append('1')
            else: printlist.append('0')
        print('\t'.join(printlist))
