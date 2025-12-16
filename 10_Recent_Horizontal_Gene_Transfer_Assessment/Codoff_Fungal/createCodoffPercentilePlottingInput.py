import os
import sys
from collections import defaultdict

larger_mapping = {'Agaricomycetes': '3. Dikarya', 'Chytridiales': '1. Chytridiomycota', 'Mucoromycota': '2. Zygomycota', 'Neocallimastigomycota': '1. Chytridiomycota', 'Pezizomycotina': '3. Dikarya', 'Pucciniomycotina': '3. Dikarya', 'Spizellomycetales and Rhizophlyctidales': '1. Chytridiomycota', 'Rhizophydiales': '1. Chytridiomycota', 'Saccharomycotina': '3. Dikarya', 'Taphrinomycotina': '3. Dikarya', 'Tremellomycetes': '3. Dikarya', 'Ustilaginomycotina': '3. Dikarya', 'Zoopagomycota': '2. Zygomycota'}

# mostly based on review by Nagy et al. 2018
cm_mapping = {'Agaricomycetes': 'Complex multicellularity', 'Chytridiales': 'Zoosporic', 'Mucoromycota': 'Some complex multicellularity', 'Neocallimastigomycota': 'Zoosporic', 'Pezizomycotina': 'Complex multicellularity', 'Pucciniomycotina': 'Some complex multicellularity', 'Spizellomycetales and Rhizophlyctidales': 'Zoosporic', 'Rhizophydiales': 'Zoosporic', 'Saccharomycotina': 'Yeast', 'Taphrinomycotina': 'Some complex multicellularity', 'Tremellomycetes': 'Some complex multicellularity', 'Ustilaginomycotina': 'Some complex multicellularity', 'Zoopagomycota': 'Other'}

ploidy_mapping = {'Agaricomycetes': 'Dikaryotic - HET enriched', 'Chytridiales': 'Diploid+ dominant', 'Mucoromycota': 'Haploid dominant', 'Neocallimastigomycota': 'Haploid dominant', 'Pezizomycotina': 'Dikaryotic - HET enriched', 'Pucciniomycotina': 'Dikaryotic', 'Spizellomycetales and Rhizophlyctidales': 'Haploid dominant', 'Rhizophydiales': 'Diploid+ dominant', 'Saccharomycotina': 'Dikaryotic', 'Taphrinomycotina': 'Dikaryotic', 'Tremellomycetes': 'Dikaryotic', 'Ustilaginomycotina': 'Dikaryotic', 'Zoopagomycota': 'Diploid+ dominant'}

reps = set([])
with open('Representatives.txt') as orf:
    for i, line in enumerate(orf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        reps.add(ls[0])

antismash_file  = 'AntiSMASH_Stats.txt'
group_mapping = {}
with open('Clade_Mapping.txt') as ocmf:
    for line in ocmf:
        line = line.strip()
        ls = line.split('\t')
        group_mapping[ls[0]] = ' '.join(ls[1].split('_'))

genome_bgc_pvals = defaultdict(list)
with open('Codoff_Results_v1.2.3.txt') as ocrf:
    for i, line in enumerate(ocrf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        g = ls[0].split('/')[1]
        bgc = ls[0].split(':Discordance')[0].split('/')[-1]
        genome_bgc_pvals[g].append(tuple([bgc, ls[-1]]))

contaminated = set(['JADGJL010000003.1.region001.gbk'])
print('\t'.join(['genome', 'bgc', 'clade', 'large_clade', 'ploidy_status', 'morphology', 'codoff_percentile']))
with open(antismash_file) as oaf:
    for i, line in enumerate(oaf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        if ls[0] not in reps:
            continue
        if not ls[0] in group_mapping: continue
        c = group_mapping[ls[0]]
        if not c in larger_mapping: continue
        lc = larger_mapping[c]
        ps = ploidy_mapping[c]
        mo = cm_mapping[c]
        for b in genome_bgc_pvals[ls[0]]:
            if b[0] not in contaminated:
                print('\t'.join([ls[1], b[0], c, lc, ps, mo, b[1]]))
