import os
import sys

clade_files = ['Pucciniomycotina.txt', 'Ustilaginomycotina.txt', 'Tremellomycetes.txt', 'Pezizomycotina.txt', 'Taphrinomycotina.txt', 'Saccharomycotina.txt', 'Rhizophydiales.txt', 'Spizellomycetales_and_Rhizophlyctidales.txt', 'Chytridiales.txt', 'Blastocladiomycota.txt', 'Neocallimastigomycota.txt', 'Zoopagomycota.txt', 'Mucoromycota.txt', 'Agaricomycetes.txt']

reps = set([])
with open('Representatives.txt') as orf:
    for line in orf:
        line = line.strip()
        ls = line.split('\t')
        reps.add(ls[0])

for f in clade_files:
    clade = f.split('.txt')[0]
    with open('Clades/' + f) as ocf:
        for line in ocf:
            line = line.strip()
            gca = '_'.join(line.split('_')[1:]).split('.')[0].replace('GCF_', 'GCA_')
            if gca in reps:
                print(gca + '\t' + clade)
