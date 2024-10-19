import os
import sys
from collections import defaultdict

clades = ['CM', 'CM', 'NCM', 'NCM']
clade_files = ['Clades/Agaricomycetes.txt', 'Clades/BGC_Enriched_Pezizomycotina.txt', 'Clades/Ascomycota_NonBGCEnrichedClade.txt', 'Clades/Basidiomycota_Yeasts_and_Other.txt']

cm_gcas = set([])
ncm_gcas = set([])
for i, f in enumerate(clade_files):
    clade = clades[i]
    with open(f) as ocf:
        for line in ocf:
            line = line.strip()
            gca = '_'.join(line.split('_')[1:]).replace('GCF_', 'GCA_').split('.')[0]
            if clade == 'CM':
                cm_gcas.add(gca)
            else:
                ncm_gcas.add(gca)

eggnog_annots = defaultdict(lambda: [0.0, (['NA']*20)])
header = []
with open('out.emapper.annotations') as oeaf:
    for line in oeaf:
        line = line.strip('\n')
        ls = line.split('\t')
        if not line.startswith('#'):
            score = float(ls[3])
            og = ls[0].split('|')[0]
            if score > eggnog_annots[og][0]:
                eggnog_annots[og] = [score, ls]
        elif line.startswith('#query'):
            header = [ls[0][1:]] + ls[1:]

with open('Significant_Hits.txt') as oshf:
    for i, line in enumerate(oshf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0:
            print('\t'.join(ls[:-3] + [ls[-1], 'CM_Conservation', 'NCM_Conservation'] + header))
        else:
            gca_with = set([x.strip() for x in ls[7].split(',')])
            cm_cons = float(len(cm_gcas.intersection(gca_with)))/len(cm_gcas)
            ncm_cons = float(len(ncm_gcas.intersection(gca_with)))/len(ncm_gcas)
            og = ls[0]
            #if (cm_cons >= 0.70 and ncm_cons <= 0.30) or (ncm_cons >= 0.70 and cm_cons <= 0.30):
            print('\t'.join(ls[:-3] + [ls[-1], str(cm_cons), str(ncm_cons)] + eggnog_annots[og][1]))
        
