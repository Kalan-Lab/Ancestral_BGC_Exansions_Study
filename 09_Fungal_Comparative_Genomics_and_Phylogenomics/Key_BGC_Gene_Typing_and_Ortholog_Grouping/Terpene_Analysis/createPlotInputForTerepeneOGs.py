import os
import sys
from collections import defaultdict
from Bio import SeqIO

pez_clade_file = 'Clades/Pezizomycotina.txt'
yea_clade_file = 'Clades/Dikarya_Yeasts.txt'
bas_clade_file = 'Clades/NonYeast_Basidi.txt'
neo_clade_file = 'Clades/Neocalli.txt'
zoo_clade_file = 'Clades/Zoopago.txt'
muc_clade_file = 'Clades/Mucoro.txt'
spi_clade_file = 'Clades/Spizello.txt'
out_clade_file = 'Clades/Other_Zoosporic_or_Basal.txt'

genus_reps_file = '../Genus_Representative_GCAs.txt'

genus_reps = set([])
with open(genus_reps_file) as ogrf:
    for line in ogrf:
        line = line.strip()
        genus_reps.add(line)


clade_lists = [pez_clade_file, bas_clade_file, yea_clade_file, muc_clade_file, zoo_clade_file, neo_clade_file, spi_clade_file, out_clade_file]
clade_names = ['1. Pezizomycotina', '2. Non-Yeast Basidimycota', '3. Yeast Dikarya', '4. Mucoromycota', '5. Zoopagomycota', '6. Neocallimastigomycota', '7. Spizellomycetales', '8. Other Zoosporic and Basal Fungi']

clade_gcas = defaultdict(set)
gca_to_clade = {}
for i, clf in enumerate(clade_lists):
    cn = clade_names[i]
    with open(clf) as ocf:
        for line in ocf:
            line = line.strip()
            gca = '_'.join(line.split('_')[-2:])
            if gca in genus_reps:
                clade_gcas[cn].add(gca)
                gca_to_clade[gca] = cn

gca_seqs = defaultdict(set)
with open('All_Key_Terpene_Proteins.faa') as oktf:
    for rec in SeqIO.parse(oktf, 'fasta'):
        rec_desc = rec.description
        gca = rec_desc.split()[1].split('|')[1]
        gca_seqs[gca].add(str(rec.seq))

select = set([])
with open('Terpene_OGs_in_atl25perGenomes.txt') as of:
    for line in of:
        line = line.strip()
        ls = line.split('\t')
        select.add(ls[0])

prot_dir = 'Proteomes/'
cdss = set([])
for f in os.listdir(prot_dir):
    gca = f.split('.faa')[0]
    with open(prot_dir + f) as opf:
        for rec in SeqIO.parse(opf, 'fasta'):
            if rec.seq in gca_seqs[gca]:
                cdss.add(rec.id)

og_gcas_bgc = defaultdict(set)
og_gcas_nonbgc = defaultdict(set)
gcas = []
with open('Orthogroups.tsv') as ogf:
    for i, line in enumerate(ogf):
        line = line.strip('\n')
        ls = line.split('\t')
        if i == 0: 
            gcas = ls[1:]
            continue
        og = ls[0]
        if og in select:
            for j, cs in enumerate(ls[1:]):
                gca = gcas[j]
                for c in cs.split(', '):
                    if c.strip() in cdss:
                        og_gcas_bgc[og].add(gca)
                    elif c.strip() != '':
                        og_gcas_nonbgc[og].add(gca)

# print('\t'.join(['ST', 'ST_Type', 'Clade', 'Prop_Clade_With', 'Pezi_Prop_Clade_With']))
clade_og_bgc = defaultdict(lambda: defaultdict(set))
for og in og_gcas_bgc:
    for gca in og_gcas_bgc[og]:
        if not gca in gca_to_clade: continue
        clade = gca_to_clade[gca]
        clade_og_bgc[og][clade].add(gca)
        
clade_og_nonbgc = defaultdict(lambda: defaultdict(set))
for og in og_gcas_nonbgc:
    for gca in og_gcas_nonbgc[og]:
        if not gca in gca_to_clade: continue
        clade = gca_to_clade[gca]
        if not gca in clade_og_bgc[og][clade]:
            clade_og_nonbgc[og][clade].add(gca)

og_name = {}
with open('OG_Annots.txt') as of:
    for line in of:
        line = line.strip()
        ls = line.split('\t')
        og_name[ls[0]] = ls[1]


st_order = {}
with open('../ST_Orders.txt') as osf:
    for line in osf:
        line = line.strip()
        ls = line.split('\t')
        st_order[ls[0]] = ls[1]

for og in og_gcas_bgc:
    for clade in clade_names:
        clade_total = len(clade_gcas[clade])
        og_bgc_context = len(clade_og_bgc[og][clade])
        og_nonbgc_context = len(clade_og_nonbgc[og][clade])
        
        og_bgc_context_prop = og_bgc_context/clade_total
        og_nonbgc_context_prop = og_nonbgc_context/clade_total

        print('\t'.join([str(x) for x in [og_name[og], 'Terpene', clade, og_bgc_context_prop, st_order[og_name[og]]]]))
        print('\t'.join([str(x) for x in [og_name[og], '0. Suspected', clade, og_nonbgc_context_prop, st_order[og_name[og]]]]))
