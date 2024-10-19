import os
import sys

larger_mapping = {'Agaricomycetes': '3. Dikarya', 'Chytridiales': '1. Chytridiomycota', 'Mucoromycota': '2. Zygomycota', 'Neocallimastigomycota': '1. Chytridiomycota', 'Pezizomycotina': '3. Dikarya', 'Pucciniomycotina': '3. Dikarya', 'Spizellomycetales and Rhizophlyctidales': '1. Chytridiomycota', 'Rhizophydiales': '1. Chytridiomycota', 'Saccharomycotina': '3. Dikarya', 'Taphrinomycotina': '3. Dikarya', 'Tremellomycetes': '3. Dikarya', 'Ustilaginomycotina': '3. Dikarya', 'Zoopagomycota': '2. Zygomycota'}

# mostly based on review by Nagy et al. 2018
cm_mapping = {'Agaricomycetes': 'Complex multicellularity', 'Chytridiales': 'Zoosporic', 'Mucoromycota': 'Some complex multicellularity', 'Neocallimastigomycota': 'Zoosporic', 'Pezizomycotina': 'Complex multicellularity', 'Pucciniomycotina': 'Some complex multicellularity', 'Spizellomycetales and Rhizophlyctidales': 'Zoosporic', 'Rhizophydiales': 'Zoosporic', 'Saccharomycotina': 'Yeast', 'Taphrinomycotina': 'Some complex multicellularity', 'Tremellomycetes': 'Some complex multicellularity', 'Ustilaginomycotina': 'Some complex multicellularity', 'Zoopagomycota': 'Other'}

ploidy_mapping = {'Agaricomycetes': 'Dikaryotic - HET enriched', 'Chytridiales': 'Diploid+ dominant', 'Mucoromycota': 'Haploid dominant', 'Neocallimastigomycota': 'Haploid dominant', 'Pezizomycotina': 'Dikaryotic - HET enriched', 'Pucciniomycotina': 'Dikaryotic', 'Spizellomycetales and Rhizophlyctidales': 'Haploid dominant', 'Rhizophydiales': 'Diploid+ dominant', 'Saccharomycotina': 'Dikaryotic', 'Taphrinomycotina': 'Dikaryotic', 'Tremellomycetes': 'Dikaryotic', 'Ustilaginomycotina': 'Dikaryotic', 'Zoopagomycota': 'Diploid+ dominant'}


result_file = 'psapas_pairwise_results/track.txt'

# apply capping
with open(result_file) as orf:
    for i, line in enumerate(orf):
        line = line.strip()
        if i == 0:
            print(line + '\tlarge_clade\tploidy_status\tmorphology\tbgcome_size')
        else:
            ls = line.split('\t')
            if ls[0] == 'Blastocladiomycota': continue
            gs = ls[2]
            c = ' '.join(ls[0].split('_')).replace('BGC Enriched Pezizomycotina', 'Pezizomycotina')
            capped_gs = min([10.0, float(gs)])
            print('\t'.join([c, ls[1], str(capped_gs)] + ls[3:] + [larger_mapping[c], ploidy_mapping[c], cm_mapping[c]]))
