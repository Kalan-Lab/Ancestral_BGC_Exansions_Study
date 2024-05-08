import os
import sys
from collections import defaultdict
from Bio import SeqIO
import json 

pez_clade_file = 'Clades/Pezizomycotina.txt'
yea_clade_file = 'Clades/Dikarya_Yeasts.txt'
bas_clade_file = 'Clades/NonYeast_Basidi.txt'
neo_clade_file = 'Clades/Neocalli.txt'
zoo_clade_file = 'Clades/Zoopago.txt'
muc_clade_file = 'Clades/Mucoro.txt'
spi_clade_file = 'Clades/Spizello.txt'
out_clade_file = 'Clades/Other_Zoosporic_or_Basal.txt'

all_pks_prots_file = 'synthaser_input_files/All_PKS_Full_Proteins.faa'
all_nrps_prots_file = 'synthaser_input_files/All_NRPS_Full_Proteins.faa'

genus_reps_file = 'Genus_Representative_GCAs.txt'

coi_set = set(['PKS|Type I|Non-reducing', 'PKS|Type I|Highly-reducing', 'PKS|Type I|Partially-reducing', 'Hybrid NRPS-PKS', 'Hybrid PKS-NRPS', 'NRPS'])
nrps_and_pks_prots_listing_file = 'NRPS_and_PKS_Proteins.txt'

pks_contam = set(['278', '279', '280', '281'])
nrps_contam = set(['144'])

redundant = set([])
redundant_names = set([])
with open(nrps_and_pks_prots_listing_file) as onaplf:
    for line in onaplf:
        line = line.strip()
        redundant.add(line)
        redundant_names.add(' '.join(line.split()[1:]))

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

clade_gca_counts = {}
for clade in clade_gcas:
    clade_gca_counts[clade] = len(clade_gcas[clade])

id_to_gca = {}
redundant_pks_ids = set([])
with open(all_pks_prots_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        if rec.id in pks_contam: continue
        gca = rec.description.split('|')[1]
        id_to_gca[rec.id] = rec.description.split('|')[1]
        if ' '.join(rec.description.split()[1:]) in redundant_names:
            redundant_pks_ids.add(rec.id)

synth_files = ['synthaser_result_files/synthaser_calls_on_pks_p1.json', 'synthaser_result_files/synthaser_calls_on_pks_p2.json', 'synthaser_result_files/synthaser_calls_on_pks_p3.json', 'synthaser_result_files/synthaser_calls_on_pks_p4.json']

classifications = set([])
clade_subtype_gcas = defaultdict(lambda: defaultdict(set))
all_sts = set([])
for sf in synth_files:
    with open(sf) as osf:
        d = json.load(osf)
        for pks in d:
            pks_id = pks['header']
            if pks_id in pks_contam: continue
            pks_gca = id_to_gca[pks_id]
            if not pks_gca in gca_to_clade: continue
            pks_clade = gca_to_clade[pks_gca]
            classified = False
            if len(pks['classification']) > 0:
                pks_class = '|'.join(pks['classification'])
                if pks_class in coi_set:
                    if pks_class == 'Hybrid PKS-NRPS':
                        pks_class = 'Hybrid NRPS-PKS'
                    classifications.add(pks_class)
                    clade_subtype_gcas[pks_clade][pks_class].add(pks_gca)
                    classified = True
            if not classified:
                class_type = 'NA'
                #if len(pks['classification']) > 0:
                #    class_type = '|'.join(pks['classification'])
                #else:
                if pks_id in redundant_pks_ids:
                    class_type = 'Suspected Hybrid-NRPS-PKS'
                else:
                    class_type = 'Suspected PKS'
                clade_subtype_gcas[pks_clade][class_type].add(pks_gca)
                classifications.add(class_type)

id_to_gca = {}
redundant_ids = set([])
with open(all_nrps_prots_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        if rec.id in nrps_contam: continue
        gca = rec.description.split('|')[1]
        if rec.description in redundant:
            redundant_ids.add(rec.id)
        id_to_gca[rec.id] = rec.description.split('|')[1]

synth_files = ['synthaser_result_files/NRPS_synthaser_p1.json', 'synthaser_result_files/NRPS_synthaser_p2.json', 'synthaser_result_files/NRPS_synthaser_p3.json']

all_sts = set([])
for sf in synth_files:
    with open(sf) as osf:
        d = json.load(osf)
        for nrps in d:
            nrps_id = nrps['header']
            if nrps_id in nrps_contam: continue
            nrps_gca = id_to_gca[nrps_id]
            if not nrps_gca in gca_to_clade: continue
            nrps_clade = gca_to_clade[nrps_gca]
            classified = False
            if len(nrps['classification']) > 0:
                nrps_class = '|'.join(nrps['classification'])
                if nrps_class in coi_set: 
                    if nrps_class == 'Hybrid PKS-NRPS':
                        nrps_class = 'Hybrid NRPS-PKS'
                    clade_subtype_gcas[nrps_clade][nrps_class].add(nrps_gca)
                    if not nrps_id in redundant_ids:
                        classified = True
                    classifications.add(nrps_class)
            if not classified and not nrps_id in redundant_ids:
                class_type = 'Suspected NRPS'
                #if len(nrps['classification']) > 0:
                #    class_type = '|'.join(nrps['classification'])
                classifications.add(class_type)
                clade_subtype_gcas[nrps_clade][class_type].add(nrps_gca)

st_order = {}
with open('ST_Orders.txt') as osf:
    for line in osf:
        line = line.strip()
        ls = line.split('\t')
        st_order[ls[0]] = ls[1]

print(st_order)

print('\t'.join(['ST', 'ST_Type', 'Clade', 'Prop_Clade_With', 'ST_Order']))
for clade in clade_names:
    for st in classifications:
        gcas_with_dt = len(clade_subtype_gcas[clade][st])
        gcas_total = 0
        
        gcas_total = clade_gca_counts[clade]
        prop_clade_with = gcas_with_dt/gcas_total


        st_name = st
        st_type = 'NA'
        if 'Suspected' in st:
            st_type = '0. Suspected'
            if st == 'Suspected NRPS':
                st_name = 'NRPS'
                gcas_with_dt = len(clade_subtype_gcas[clade][st].difference(clade_subtype_gcas[clade]['NRPS']))
            elif st == 'Suspected Hybrid-NRPS-PKS':
                st_name = 'Hybrid NRPS-PKS'
                gcas_with_dt = len(clade_subtype_gcas[clade][st].difference(clade_subtype_gcas[clade]['Hybrid NRPS-PKS']))
            else:
                st_name = 'Only Suspected PKS'
                gcas_with_dt = len(clade_subtype_gcas[clade][st].difference(clade_subtype_gcas[clade]['PKS|Type I|Non-reducing']).difference(clade_subtype_gcas[clade]['PKS|Type I|Highly-reducing']).difference(clade_subtype_gcas[clade]['PKS|Type I|Partially-reducing']))

        elif 'NRPS' in st and 'PKS' in st:
            st_type = 'Hybrid'
        elif 'NRPS' in st:
            st_type = 'NRPS'
        elif 'PKS' in st:
            st_type = 'PKS'
        
        prop_clade_with = gcas_with_dt/gcas_total
        order = st_order[st_name]
        print('\t'.join([str(x) for x in [st_name, st_type, clade, prop_clade_with, order]]))
