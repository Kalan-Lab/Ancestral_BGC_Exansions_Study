import os
import sys
import numpy

def p_adjust_bh(p):
    """
    Using implementation from https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python/33532498#33532498
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    """
    p = numpy.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / numpy.arange(len(p), 0, -1)
    q = numpy.minimum(1, numpy.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

permutation_file = sys.argv[1]

data = []
pvals = []
with open(permutation_file) as opf:
    for i, line in enumerate(opf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0:
            print('\t'.join(['adj_pval'] + ls))
            continue
        pval = float(ls[1])
        pvals.append(pval)
        data.append(ls)

adj_pvals = p_adjust_bh(pvals)

for i, aj in enumerate(adj_pvals):
    if aj < 0.05:
        print(str(aj) + '\t' + '\t'.join(data[i]))
