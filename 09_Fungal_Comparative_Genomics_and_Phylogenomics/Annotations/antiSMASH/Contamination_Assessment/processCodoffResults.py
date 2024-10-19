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

alpha = 0.1

data = []
pvals = []
input_file = 'Codoff_Results.txt'
with open(input_file) as oif:
    for i, line in enumerate(oif):
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('/')[1]
        pval = float(ls[1])
        bgc = ls[0].split('/')[2].split(':Empirical P-value')[0]
        data.append([gca, bgc, pval])
        pvals.append(pval)

adj_pvals = p_adjust_bh(pvals)
for i, ls in enumerate(data):
    adj_pval = adj_pvals[i]

    if adj_pval < alpha:
        print('\t'.join([str(x) for x in (ls + [adj_pval])]))

