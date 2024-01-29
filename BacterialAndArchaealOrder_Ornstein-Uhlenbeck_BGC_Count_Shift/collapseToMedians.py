import os
import sys
from collections import defaultdict
import statistics

"""
0       Assembly_ID
1       Order
2       Genus
3       BGC_Count
4       BGCome_Size
5       Genome_Size
6       BGC_Genome_Prop
7       Metallophore_Size
8       NRPS_Size
9       PKS_Size
10      Terpene_Size
11      Distinct_KOfams
12      Distinct_CAZy
"""

overview_file = 'Overview_GCA_Level.txt'
treeorder_file = 'Tree_Order.txt'

gs_list = defaultdict(list)
bs_list = defaultdict(list)
ts_list = defaultdict(list)
kf_list = defaultdict(list)
cf_list = defaultdict(list)
ns_list = defaultdict(list)
ps_list = defaultdict(list)
orders = set([])
with open(overview_file) as ovf:
    for i, line in enumerate(ovf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        order = ls[1]
        gs = int(ls[5])
        bs = float(ls[5])*float(ls[6])
        ts = float(ls[10])
        ns = float(ls[8])
        ps = float(ls[9])
        kf = float(ls[11])
        cf = float(ls[12])
        gs_list[order].append(gs)
        bs_list[order].append(bs)
        ts_list[order].append(ts)
        ns_list[order].append(ns)
        ps_list[order].append(ps)
        kf_list[order].append(kf)
        cf_list[order].append(cf)
        orders.add(order)

data = {}
print('\t'.join(['order', 'median_genome_size', 'median_bgc_size', 'median_terpene_count', 'median_nrps_count', 'median_pks_count', 'median_oxy_count', 'median_cazy_count']))
for o in orders:
    mgs = statistics.median(gs_list[o])
    mbs = statistics.median(bs_list[o])
    mts = statistics.median(ts_list[o])
    mns = statistics.median(ns_list[o])
    mps = statistics.median(ps_list[o])
    mkf = statistics.median(kf_list[o])
    mcf = statistics.median(cf_list[o])
    row = '\t'.join([str(x) for x in [o, mgs, mbs, mts, mns, mps, mkf, mcf]])
    data[o] = row

with open(treeorder_file) as otf:
    for line in otf:
        line = line.strip()
        ls = line.split('\t')
        print(data[line])
