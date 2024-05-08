import os
import sys
from collections import defaultdict
import statistics
import numpy as np

overview_gca_file = 'Overview_GCA_Level.txt'
itol_dir = 'iTol_Tracks/'
factors = []
orders = set([])
factor_order_values = defaultdict(lambda: defaultdict(list))
with open(overview_gca_file) as ogf:
    for i, line in enumerate(ogf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: 
            factors = ls[3:]
        else:
            order = ls[1]
            orders.add(order)
            for j, val in enumerate(ls[3:]):
                factor = factors[j]
                val = float(val)
                factor_order_values[factor][order].append(val)

for f in factors:
    outf = open(itol_dir + f + '.iTol.txt', 'w')
    outf.write('DATASET_BOXPLOT\n')
    outf.write('SEPARATOR TAB\n')
    outf.write('COLOR\t#000000\n')
    outf.write('DATASET_LABEL\t' + f + '\n')  
    outf.write('DATA\n')
    for o in orders:
        tbc = factor_order_values[f][o]
        low95 = np.percentile(tbc, 5)
        quart25 = np.percentile(tbc, 25)
        med = np.percentile(tbc, 50)
        quart75 = np.percentile(tbc, 75)
        high95 = np.percentile(tbc, 95)
        outliers = []
        for bs in tbc:
            if bs < low95 or bs > high95:
                outliers.append(bs)
        outf.write(o + '\t' + '\t'.join([str(x) for x in ([low95, quart25, med, quart75, high95] + outliers)]) + '\n')
    outf.close()
