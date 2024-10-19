import os
import sys
from collections import defaultdict

pez_clade_file = 'Clades/Pezizomycotina.txt'
asc_clade_file = 'Clades/Ascomycota_Yeasts_and_Other.txt'
bas_clade_file = 'Clades/Basidiomycota_Yeasts_and_Other.txt'
aga_clade_file = 'Clades/Agaricomycetes.txt'
neo_clade_file = 'Clades/Neocallimastigomycota.txt'
zoo_clade_file = 'Clades/Zoopagomycota.txt'
muc_clade_file = 'Clades/Mucoromycota.txt'
chy_clade_file = 'Clades/Chytrid_and_ChytridLike.txt'
bla_clade_file = 'Clades/Blastocladiomycota.txt'
out_clade_file = 'Clades/Basal_Fungi.txt'

genus_reps_file = 'Representatives.txt'

genus_reps = set([])
with open(genus_reps_file) as ogrf:
    for line in ogrf:
        line = line.strip()
        genus_reps.add(line.split('\t')[0])

clade_lists = [pez_clade_file, asc_clade_file, aga_clade_file, bas_clade_file, muc_clade_file, zoo_clade_file, neo_clade_file, chy_clade_file, bla_clade_file, out_clade_file]
clade_names = ['01. Pezizomycotina', '02. Yeast and Other Ascomycota', '03. Agaricomycetes', '04. Yeast and Other Basidiomycota', '05. Mucoromycota', '06. Zoopagomycota', '07. Neocallimastigomycota', '08. Chytrid and Other Chytrid-Like', '09. Blastocladiomycota', '10. Basal Fungi']

clade_gcas = defaultdict(set)
gca_to_clade = {}
for i, clf in enumerate(clade_lists):
    cn = clade_names[i]
    with open(clf) as ocf:
        for line in ocf:
            line = line.strip()
            gca = '_'.join(line.split('_')[1:]).replace('GCF_', 'GCA_').split('.')[0]
            if gca in genus_reps:
                clade_gcas[cn].add(gca)
                gca_to_clade[gca] = cn

clade_gca_counts = {}
for clade in clade_gcas:
    clade_gca_counts[clade] = len(clade_gcas[clade])

dom_gcas = defaultdict(set)
with open('Results.txt') as orf:
    for line in orf:
        line = line.strip()
        ls = line.split('\t')

        dom = ls[0]
        gca = ls[2].split('|')[0]
        if gca in genus_reps:
            dom_gcas[dom].add(gca)

order = {}
with open('ST_Orders.txt') as osof:
    for line in osof:
        line = line.strip()
        ls = line.split('\t')
        order[ls[0]] = ls[1]

print('\t'.join(['ST', 'ST_Type', 'Clade', 'Prop_Clade_With', 'ST_Order']))
for dom in dom_gcas:
    for clade in clade_gcas:
        dom_gcas_in_clade = clade_gcas[clade].intersection(dom_gcas[dom])
        prop_clade_with = len(dom_gcas_in_clade)/float(len(clade_gcas[clade]))
        print('\t'.join([dom, 'Terpene', clade, str(prop_clade_with), order[dom]]))
