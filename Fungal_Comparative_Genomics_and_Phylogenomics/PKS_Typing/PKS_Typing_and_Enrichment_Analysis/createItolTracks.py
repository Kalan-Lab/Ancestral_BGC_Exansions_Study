import os
import sys
from collections import defaultdict
from ete3 import Tree

resdir = 'Results/'
tree_file = '../KS_Proteins.msa.trimal_strict.fasttree_lg.tre'
itol_dir = 'iTol_Tracks/'

red_sup = {}
nonred_sup = {}
pks_type = {}
nrps_part = {}
for b in os.listdir(resdir):
    pks_file = resdir + b + '/pks_type.txt'
    with open(pks_file) as opf:
        for line in opf:
            line = line.strip('\n')
            ls = line.split('\t')
            pks_type[b] = ls[0]
            if ls[0] == '':
                pks_type[b] = 'NA'
            nrps_part[b] = ls[3]
            red_sup[b] = ls[1]
            nonred_sup[b] = ls[2]

pks_type_colors = {'HR (highly-reducing)': '#eb9534', 'NR (non-reducing)': '#ede7b7', 'NA': '#bababa'}
part_nrps_colors = {'False': '#ffffff', 'True': '#a64d79'}

outf1 = open(itol_dir + 'pks_status.itol.txt', 'w')
outf2 = open(itol_dir + 'part_nrps.itol.txt', 'w')
outf3 = open(itol_dir + 'reducing_support.itol.txt', 'w')
outf4 = open(itol_dir + 'nonreducing_support.itol.txt', 'w')

outf1.write('DATASET_COLORSTRIP\n')
outf1.write('SEPARATOR TAB\n')
outf1.write('DATASET_LABEL\tpks_status\n')
outf1.write('COLOR\t#000000\n')
outf1.write('DATA\n')

outf2.write('DATASET_COLORSTRIP\n')
outf2.write('SEPARATOR TAB\n')
outf2.write('DATASET_LABEL\tpart_nrps\n')
outf2.write('COLOR\t#000000\n')
outf2.write('DATA\n')

outf3.write('DATASET_SIMPLEBAR\n')
outf3.write('SEPARATOR TAB\n')
outf3.write('DATASET_LABEL\treducing_support\n')
outf3.write('COLOR\t#000000\n')
outf3.write('DATA\n')

outf4.write('DATASET_SIMPLEBAR\n')
outf4.write('SEPARATOR TAB\n')
outf4.write('DATASET_LABEL\tnon-reducing_support\n')
outf4.write('COLOR\t#000000\n')
outf4.write('DATA\n')

t = Tree(tree_file) 
for n in t.traverse('postorder'):
    if n.is_leaf():
        name = n.name
        key = '_bgc_'.join(name.split('|')[:2])
        pt = pks_type[key]
        pn = nrps_part[key]
        rs = red_sup[key]
        ns = nonred_sup[key]
        outf1.write(name + '\t' + pks_type_colors[pt] + '\t' + pt + '\n')
        outf2.write(name + '\t' + part_nrps_colors[pn] + '\t' + pn + '\n')
        outf3.write(name + '\t' + rs + '\n')
        outf4.write(name + '\t' + ns + '\n')
outf1.close()
outf2.close()
outf3.close()
outf4.close()
