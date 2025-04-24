import os
import sys
from collections import defaultdict
import statistics

"""
0	Assembly_ID
1	Order
2	Genus
3	BGC_Count
4	BGCome_Size
5	Genome_Size
6	BGC_Genome_Prop
7	Metallophore_Size
8	NRPS_Size
9	PKS_Size
10	NRPS_or_PKS_Size
11	Relaxed_BGCome_Size
12	Relaxed_BGC_Count
13	Distinct_KOfams
14	Distinct_CAZy
15	Total_CAZy
"""

overview_file = 'Overview_GCA_Level.txt'
treeorder_file = 'Tree_Order.txt'

gs_list = defaultdict(list)
bc_list = defaultdict(list)
bs_list = defaultdict(list)
rbs_list = defaultdict(list)
rbc_list = defaultdict(list)
kf_list = defaultdict(list)
cf_list = defaultdict(list)
ns_list = defaultdict(list)
ps_list = defaultdict(list)
nps_list = defaultdict(list)
tcf_list = defaultdict(list)
orders = set([])
with open(overview_file) as ovf:
    for i, line in enumerate(ovf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        order = ls[1]
        gs = int(ls[5])
        bs = float(ls[4])
        bc = int(ls[3])
        rbs = float(ls[11])
        rbc = float(ls[12])
        ns = float(ls[8])
        ps = float(ls[9])
        nps = float(ls[10])
        kf = float(ls[13])
        cf = float(ls[14])
        tcf = float(ls[15])
        gs_list[order].append(gs)
        bs_list[order].append(bs)
        bc_list[order].append(bc)
        rbs_list[order].append(rbs)
        rbc_list[order].append(rbc)
        nps_list[order].append(nps)
        ns_list[order].append(ns)
        ps_list[order].append(ps)
        kf_list[order].append(kf)
        cf_list[order].append(cf)
        tcf_list[order].append(tcf)
        orders.add(order)

data = {}
print('\t'.join(['order', 'median_genome_size', 'median_bgcome_size', 'median_bgc_count', 'median_bgcome_size_relaxed', 'median_bgc_count_relaxed', 'median_npbgcome_size', 'median_nrps_count', 'median_pks_count', 'median_oxy_count', 'median_distinct_cazy_count', 'median_total_cazy_count']))
for o in orders:
    mgs = statistics.median(gs_list[o])
    mbs = statistics.median(bs_list[o])
    mrbs = statistics.median(rbs_list[o])
    mrbc = statistics.median(rbc_list[o])
    mns = statistics.median(ns_list[o])
    mps = statistics.median(ps_list[o])
    mnps = statistics.median(nps_list[o])
    mkf = statistics.median(kf_list[o])
    mcf = statistics.median(cf_list[o])
    mtcf = statistics.median(tcf_list[o])
    mbc = statistics.median(bc_list[o])
    row = '\t'.join([str(x) for x in [o, mgs, mbs, mbc, mrbs, mrbc, mnps, mns, mps, mkf, mcf, mtcf]])
    data[o] = row

with open(treeorder_file) as otf:
    for line in otf:
        line = line.strip()
        ls = line.split('\t')
        print(data[line])
