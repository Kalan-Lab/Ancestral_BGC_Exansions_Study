import os
import sys
from collections import defaultdict
from scipy import stats

bgc_rich = set([])
all_fungi = set([x.strip() for x in open('All_Fungi.txt').readlines()])
for cf in ['Agaricomycetes.txt', 'BGC_Enriched_Pezizomycotina.txt']:
    with open(cf) as ocf:
        for line in ocf:
            line = line.strip()
            bgc_rich.add(line)

genus_reps = set([])
with open('Representatives.txt') as orf:
    for line in orf:
        line = line.strip()
        ls = line.split('\t')
        genus_reps.add(ls[0])

prot_count = {}
with open('Phylogenetic_Regression_Input.txt') as of:
    for i, line in enumerate(of):
        if i == 0: continue
        line = line.strip('\n')
        ls = line.split('\t')
        prot_count[ls[0]] = float(ls[-1])

print('group\tgenome\thi_prop')
bgcrich = []
other = []
with open('HET_Hits.iTol.txt') as ohf:
    for i, line in enumerate(ohf):
        if i > 5:
            line = line.strip()
            ls = line.split('\t')
            gca = '_'.join(ls[0].split('_')[1:]).replace('GCF_', 'GCA_').split('.')[0]
            if not gca in genus_reps: continue
            group = None
            tp = prot_count[ls[0]]
            hp = int(ls[1]) + int(ls[2])
            hprop = hp/tp
            if ls[0] in bgc_rich:
                group = 'Agaricomycetes or BGC-rich Pezizomycotina'
                bgcrich.append(hprop)
            elif ls[0] in all_fungi:
                group = 'Other fungi'
                other.append(hprop)
            else: continue
            print(group + '\t' + ls[0] + '\t' + str(hprop))
stat, pval = stats.ranksums(bgcrich, other, alternative='greater')
sys.stderr.write(str(pval) + '\n')
