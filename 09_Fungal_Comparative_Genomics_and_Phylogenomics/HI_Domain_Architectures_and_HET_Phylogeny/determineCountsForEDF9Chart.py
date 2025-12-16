import os
import sys
from collections import defaultdict
import statistics
import numpy as np

agarico = set(['_'.join(x.strip().split('_')[-2:]).split('.')[0].replace('GCF_', 'GCA_') for x in open('HET_Domain_Phylo/Clade_Info/Agaricomycetes.txt').readlines()])
bgc_rich_pezi = set(['_'.join(x.strip().split('_')[-2:]).split('.')[0].replace('GCF_', 'GCA_') for x in open('HET_Domain_Phylo/Clade_Info/BGC_Enriched_Pezizomycotina.txt')])
all_fungi = set(['_'.join(x.strip().split('_')[-2:]).split('.')[0].replace('GCF_', 'GCA_') for x in open('HET_Domain_Phylo/Clade_Info/All_Fungi.txt').readlines()])

genus_reps = set([])
with open('Representatives.txt') as orf:
    for line in orf:
        line = line.strip()
        ls = line.split('\t')
        genus_reps.add(ls[0])

agarico = agarico.intersection(genus_reps)
bgc_rich_pezi = bgc_rich_pezi.intersection(genus_reps)
all_fungi = all_fungi.intersection(genus_reps)

prot_count = {}
with open('Phylogenetic_Regression_Input.txt') as of:
    for i, line in enumerate(of):
        if i == 0: continue
        line = line.strip('\n')
        ls = line.split('\t')
        prot_count[ls[0]] = float(ls[-1])

dom_architect = ['HET', 'Het-C', 'ATP-cone|Ribonuc_red_lgN|Ribonuc_red_lgC', 'HMG_box', 'MAT1-1-2', 'Patatin|NB-ARC', 'HeLo|HET-s_218-289', 'HET|NACHT|WD40', 'SET|Rubis-subs-bind|zf-MYND', 'Peptidase_S8', 'GLTP', 'RVT_2', 'SERF-like_N', 'CYSTM', 'Ccdc124', 'Fructosamin_kin', 'Patatin|TPR_10']
das = set(dom_architect)

print('dom_architecture\tbgc_rich_pezi\tagarico\tother')

group_gca_dom_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
with open('Domain_Architectures_of_Proteins_with_OneOrMoreHIAssociatedDoms.txt') as odf:
    for line in odf:
        line = line.strip()
        ls = line.split('\t')
        dom_arc = ls[1]
        if dom_arc in das:
            gca = ls[0].split('|')[0]
            if gca in bgc_rich_pezi:
                group_gca_dom_counts['brp'][gca][dom_arc] += 1
            elif gca in agarico:
                group_gca_dom_counts['aga'][gca][dom_arc] += 1
            else:
                group_gca_dom_counts['other'][gca][dom_arc] += 1

for da in dom_architect:
    brp_counts = []
    aga_counts = []
    other_counts = []

    for g in bgc_rich_pezi:
        brp_counts.append(group_gca_dom_counts['brp'][g][da])
    for g in agarico:
        aga_counts.append(group_gca_dom_counts['aga'][g][da])
    for g in all_fungi:
        if g in bgc_rich_pezi or g in agarico: continue
        other_counts.append(group_gca_dom_counts['other'][g][da])

    brp_mean = str(statistics.mean(brp_counts))
    brp_lq = str(np.percentile(sorted(brp_counts), 25))
    brp_uq = str(np.percentile(sorted(brp_counts), 75))

    aga_mean = str(statistics.mean(aga_counts))
    aga_lq = str(np.percentile(sorted(aga_counts), 25))
    aga_uq = str(np.percentile(sorted(aga_counts), 75))

    other_mean = str(statistics.mean(other_counts))
    other_lq = str(np.percentile(sorted(other_counts), 25))
    other_uq = str(np.percentile(sorted(other_counts), 75))


    print(da + '\t' + brp_mean + ' (' + brp_lq + '-' + brp_uq + ')\t' + aga_mean + ' (' + aga_lq + '-' + aga_uq + ')\t' + other_mean + ' (' + other_lq  + '-' + other_uq + ')')
