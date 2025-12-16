import os
import sys
from collections import defaultdict

bgc_listing_file = 'BGC_Information.txt'
genomad_dir = 'geNomad_Results/'
as_stats_file = 'AntiSMASH_Stats.Updated.txt'

sample_mge_coords = defaultdict(lambda: defaultdict(set))
for j, s in enumerate(os.listdir(genomad_dir)):
    if j%1000 == 0:
        sys.stderr.write(str(j) + '\n')
    sample = s.split('.')[0]
    vir_sum_file = genomad_dir + s + '/' + s + '_summary/' + s + '_virus_summary.tsv'
    pla_sum_file = genomad_dir + s + '/' + s + '_summary/' + s + '_plasmid_summary.tsv'

    if os.path.isfile(vir_sum_file):
        with open(vir_sum_file) as ovsf:
            for i, line in enumerate(ovsf):
                if i == 0: continue
                line = line.strip()
                ls = line.split('\t')
                scaff = ls[0].split('|')[0]
                coords = ls[3]
                if coords == 'NA':
                    length =int(ls[1])
                    sample_mge_coords[sample][scaff].add(tuple([1, length+1]))
                else:
                    cs = [int(x) for x in coords.split('-')]
                    sample_mge_coords[sample][scaff].add(tuple([cs[0], cs[1]+1]))

    if os.path.isfile(pla_sum_file):
        with open(pla_sum_file) as opsf:
            for i, line in enumerate(opsf):
                if i == 0: continue
                line = line.strip()
                ls = line.split('\t')
                scaff = ls[0].split('|')[0]
                length =int(ls[1])
                sample_mge_coords[sample][scaff].add(tuple([1, length+1]))

gca_to_phylum = {}
with open(as_stats_file) as oasf:
    for i, line in enumerate(oasf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca_to_phylum[ls[0].split('.')[0]] = ls[1]

with open(bgc_listing_file) as oblf:
    for line in oblf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('.')[0]
        phylum = gca_to_phylum[gca]
        overlap_with_mge = False
        scaff = ls[2]
        start = int(ls[3])
        end = int(ls[4])
        bgc_coords = set(range(start+1, end+1))
        mge_coords = set([])
        for mc in sample_mge_coords[gca][scaff]:
            mc = list(mc)
            mge_coords = mge_coords.union(set(range(mc[0], mc[1])))
        if len(bgc_coords.intersection(mge_coords)) > 0:
            overlap_with_mge = True
        bgc = ls[1]
        print(gca + '\t' + bgc + '\t' + phylum + '\t' + str(overlap_with_mge))
