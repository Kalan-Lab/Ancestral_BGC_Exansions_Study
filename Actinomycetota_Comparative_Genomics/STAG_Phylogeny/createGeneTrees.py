import os
import sys

prot_dir = os.path.abspath(sys.argv[1])+'/'
tree_dir = os.path.abspath(sys.argv[2])+'/'

for f in os.listdir(prot_dir):
    gn = f.split('.fa')[0]
    gn_dir = tree_dir + gn + '/'
    os.system('mkdir ' + gn_dir)
    print('iqtree2 -s ' + prot_dir + f + ' -mset LG,WAG --threads-max 5 -B 1000  -pre ' + gn_dir + gn)
