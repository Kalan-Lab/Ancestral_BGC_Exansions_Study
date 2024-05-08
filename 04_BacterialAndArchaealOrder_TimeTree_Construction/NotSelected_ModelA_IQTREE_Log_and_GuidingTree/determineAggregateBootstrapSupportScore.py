import os
import sys
from ete3 import Tree
import statistics

bss = []
t = Tree(sys.argv[1])
for n in t.traverse():
    if not n.is_leaf():
        bss.append(n.support)

print('Sum:\t' + str(sum(bss)))
print('Median:\t' + str(statistics.median(bss)))
print('Num nodes with support:\t' + str(bss))
