import os
import sys
from ete3 import Tree

t = Tree(sys.argv[1])

print('DATASET_COLORSTRIP')
print('SEPARATOR TAB')
print('COLOR\t#000000')
print('DATASET_LABEL\tMIBIG_or_Basal')
print('DATA')
for n in t.traverse('postorder'):
    name = n.name
    if n.is_leaf():
        if 'MIBiG' in name:
            print(name + '\t#072375\tMIBiG')
        else:
            print(name + '\t#818285\tBasal')
