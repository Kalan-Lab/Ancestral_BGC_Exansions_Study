import os
import sys
from collections import defaultdict
from Bio import SeqIO
import json 

pez_clade_file = 'Taxonomic_Clades/Pezizomycotina.txt'
mul_clade_file = 'Taxonomic_Clades/Multicellular_Clade.txt'
asc_clade_file = 'Taxonomic_Clades/Ascomycota.txt'
zyg_clade_file = 'Taxonomic_Clades/Zygomycota.txt'
bas_clade_file = 'Taxonomic_Clades/Basidiomycota.txt'
all_fungi_file = 'Taxonomic_Clades/All_Genomes.txt'
all_pks_prots_file = '../All_PKS_Full_Proteins.faa'
all_nrps_prots_file = '../All_NRPS_Full_Proteins.faa'
genus_reps_file = 'Genus_Representative_GCAs.txt'
coi_set = set(['PKS|Type I|Non-reducing', 'PKS|Type I|Highly-reducing', 'PKS|Type I|Partially-reducing', 'Hybrid NRPS-PKS', 'Hybrid PKS-NRPS', 'NRPS'])
#coi_set = set(['Fatty acid synthase', 'Fatty acid synthase|alpha subunit'])
nrps_and_pks_prots_listing_file = 'NRPS_and_PKS_Proteins.txt'

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

asc_clade = set([])
mul_clade = set([])
pez_clade = set([])
all_fungi = set([])
bas_clade = set([])
zyg_clade = set([])

with open(all_fungi_file) as oaff:
    for line in oaff:
        line = line.strip()
        all_fungi.add(line)

with open(pez_clade_file) as obef:
    for line in obef:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        pez_clade.add(gca)

with open(mul_clade_file) as odcf:
    for line in odcf:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        mul_clade.add(gca)

with open(bas_clade_file) as obcf:
    for line in obcf:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        bas_clade.add(gca)

with open(zyg_clade_file) as ozcf:
    for line in ozcf:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        zyg_clade.add(gca)

with open(asc_clade_file) as oacf:
    for line in oacf:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        asc_clade.add(gca)

other_fungi = all_fungi.intersection(genus_reps).difference(mul_clade)
mul_clade = mul_clade.intersection(genus_reps)
pez_clade = pez_clade.intersection(genus_reps)
asc_clade = asc_clade.intersection(genus_reps).difference(pez_clade)
bas_clade = bas_clade.intersection(genus_reps)
zyg_clade = zyg_clade.intersection(genus_reps)

id_to_gca = {}
redundant_pks_ids = set([])
clade_counts = defaultdict(int)
clade_with_putative_pks_counts = defaultdict(set)
with open(all_pks_prots_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        gca = rec.description.split('|')[1]
        clade = None
        if gca in pez_clade:
            clade = '1. Pezizomycotina'
        elif gca in asc_clade:
            clade = '2. Other Ascomycota'
        elif gca in bas_clade:
            clade = '3. Basidiomycota'
        elif gca in zyg_clade:
            clade = '4. Zygomycota'
        elif gca in other_fungi:
            clade = '5. Zoosporic Fungi'
        clade_counts[clade] += 1
        id_to_gca[rec.id] = rec.description.split('|')[1]
        clade_with_putative_pks_counts[clade].add(gca)
        if ' '.join(rec.description.split()[1:]) in redundant_names:
            redundant_pks_ids.add(rec.id)

synth_files = ['synthaser_result_files/synthaser_calls_on_pks_p1.json', 'synthaser_result_files/synthaser_calls_on_pks_p2.json', 'synthaser_result_files/synthaser_calls_on_pks_p3.json', 'synthaser_result_files/synthaser_calls_on_pks_p4.json']

clade_subtype_gcas = defaultdict(lambda: defaultdict(set))
clade_classified = defaultdict(set)
unclassified_clade_counts = defaultdict(lambda: defaultdict(int))
all_sts = set([])
for sf in synth_files:
    with open(sf) as osf:
        d = json.load(osf)
        for pks in d:
            pks_id = pks['header']
            pks_gca = id_to_gca[pks_id]
            pks_clade = None
            if pks_gca in pez_clade:
                pks_clade = '1. Pezizomycotina'
            elif pks_gca in asc_clade:
                pks_clade = '2. Other Ascomycota'
            elif pks_gca in bas_clade:
                pks_clade = '3. Basidiomycota'
            elif pks_gca in zyg_clade:
                pks_clade = '4. Zygomycota'
            elif pks_gca in other_fungi:
                pks_clade = '5. Zoosporic Fungi'
            classified = False
            if pks_clade != None and len(pks['classification']) > 0:
                pks_class = '|'.join(pks['classification'])
                if pks_class in coi_set:
                    if pks_class == 'Hybrid PKS-NRPS':
                        pks_class = 'Hybrid NRPS-PKS'
                    clade_subtype_gcas[pks_clade][pks_class].add(pks_gca)
                    clade_classified[pks_clade].add('pks|' + pks_id)
                    classified = True
            if not classified and pks_clade != None:
                class_type = 'NA'
                if len(pks['classification']) > 0:
                    class_type = '|'.join(pks['classification'])
                else:
                    if pks_id in redundant_pks_ids:
                        class_type = 'NA|Suspected PKS & NRPS'
                    else:
                        class_type = 'NA|Suspected PKS'
                unclassified_clade_counts[pks_clade][class_type] += 1

id_to_gca = {}
redundant_ids = set([])
clade_with_putative_nrps_counts = defaultdict(set)
with open(all_nrps_prots_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        gca = rec.description.split('|')[1]
        clade = None
        if gca in pez_clade:
            clade = '1. Pezizomycotina'
        elif gca in asc_clade:
            clade = '2. Other Ascomycota'
        elif gca in bas_clade:
            clade = '3. Basidiomycota'
        elif gca in zyg_clade:
            clade = '4. Zygomycota'
        elif gca in other_fungi:
            clade = '5. Zoosporic Fungi'
        if not rec.description in redundant:
            clade_counts[clade] += 1    
        else:
            redundant_ids.add(rec.id)
        id_to_gca[rec.id] = rec.description.split('|')[1]
        clade_with_putative_nrps_counts[clade].add(gca)

synth_files = ['synthaser_result_files/NRPS_synthaser_p1.json', 'synthaser_result_files/NRPS_synthaser_p2.json', 'synthaser_result_files/NRPS_synthaser_p3.json']

all_sts = set([])
for sf in synth_files:
    with open(sf) as osf:
        d = json.load(osf)
        for nrps in d:
            nrps_id = nrps['header']
            nrps_gca = id_to_gca[nrps_id]
            nrps_clade = None
            if nrps_gca in pez_clade:
                nrps_clade = '1. Pezizomycotina'
            elif nrps_gca in asc_clade:
                nrps_clade = '2. Other Ascomycota'
            elif nrps_gca in bas_clade:
                nrps_clade = '3. Basidiomycota'
            elif nrps_gca in zyg_clade:
            	nrps_clade = '4. Zygomycota'
            elif nrps_gca in other_fungi:
                nrps_clade = '5. Zoosporic Fungi'
            
            classified = False
            if nrps_clade != None and len(nrps['classification']) > 0:
                nrps_class = '|'.join(nrps['classification'])
                if nrps_class in coi_set: 
                    if nrps_class == 'Hybrid PKS-NRPS':
                        nrps_class = 'Hybrid NRPS-PKS'

                    clade_subtype_gcas[nrps_clade][nrps_class].add(nrps_gca)
                    if not nrps_id in redundant_ids:
                        clade_classified[nrps_clade].add('nrps|' + nrps_id)
                        classified = True

            if not classified and not nrps_id in redundant_ids and nrps_clade != None:
                class_type = 'NA|Suspected NRPS'
                if len(nrps['classification']) > 0:
                    class_type = '|'.join(nrps['classification'])
                unclassified_clade_counts[nrps_clade][class_type] += 1

# Some of this is commented out because it was used for Supplementary figure generation
"""
print('clade\tcategory\tcount')
for clade in ['1. Pezizomycotina', '2. Other Ascomycota',  '3. Basidiomycota', '4. Zygomycota', '5. Zoosporic Fungi']:
    for class_type in unclassified_clade_counts[clade]:
        print('\t'.join([clade, class_type, str(unclassified_clade_counts[clade][class_type])]))
"""

"""
print('clade\ttype\tfrequency')
for clade in ['1. Pezizomycotina', '2. Other Ascomycota',  '3. Basidiomycota', '4. Zygomycota', '5. Zoosporic Fungi']:
    gcas_total = 0
    if clade == '1. Pezizomycotina':
        gcas_total = len(pez_clade)
    if clade == '2. Other Ascomycota':
        gcas_total = len(asc_clade)
    elif clade == '3. Basidiomycota':
        gcas_total = len(bas_clade)
    elif clade == '4. Zygomycota':
        gcas_total = len(zyg_clade)
    elif clade == '5. Zoosporic Fungi':
        gcas_total = len(other_fungi)
    gca_with_pks_like = len(clade_with_putative_pks_counts[clade])
    gca_with_nrps_like = len(clade_with_putative_nrps_counts[clade])
    print(clade + '\tPKS(-like)\t' + str(gca_with_pks_like/gcas_total))
    print(clade + '\tNRPS(-like)\t' + str(gca_with_nrps_like/gcas_total))
"""

"""
print('clade\tcategory\tcount')
for clade in ['1. Pezizomycotina', '2. Other Ascomycota',  '3. Basidiomycota', '4. Zygomycota', '5. Zoosporic Fungi']:
    clade_count = clade_counts[clade]
    classified = len(clade_classified[clade])
    print(clade + '\tclassified\t' + str(classified))
    print(clade + '\tunclassified\t' + str(clade_count-classified))
"""

"""
print('\t'.join(['clade', 'prop_clade_with_fas']))
for clade in ['1. Pezizomycotina', '2. Other Ascomycota',  '3. Basidiomycota', '4. Zygomycota', '5. Zoosporic Fungi']:
    all_gcas_with_fas = set([])
    for st in sorted(coi_set):
        all_gcas_with_fas = all_gcas_with_fas.union(clade_subtype_gcas[clade][st])

    num_gcas_with_fas = len(all_gcas_with_fas)
    gcas_total = 0

    if clade == '1. Pezizomycotina':
        gcas_total = len(pez_clade)
    elif clade == '2. Other Ascomycota':
        gcas_total = len(asc_clade)
    elif clade == '3. Basidiomycota':
        gcas_total = len(bas_clade)
    elif clade == '4. Zygomycota':
        gcas_total = len(zyg_clade)
    elif clade == '5. Zoosporic Fungi':
        gcas_total = len(other_fungi)

    prop_clade_with = num_gcas_with_fas/gcas_total
    st_type = 'FAS'
    print('\t'.join([str(x) for x in [clade, prop_clade_with]]))
"""


for clade in ['1. Pezizomycotina', '2. Other Ascomycota',  '3. Basidiomycota', '4. Zygomycota', '5. Zoosporic Fungi']:
    for st in sorted(coi_set):
        if st == 'Hybrid PKS-NRPS': continue
        gcas_with_dt = len(clade_subtype_gcas[clade][st])
        gcas_total = 0
        if clade == '1. Pezizomycotina':
            gcas_total = len(pez_clade)
        elif clade == '2. Other Ascomycota':
            gcas_total = len(asc_clade)
        elif clade == '3. Basidiomycota':
            gcas_total = len(bas_clade)
        elif clade == '4. Zygomycota':
            gcas_total = len(zyg_clade)
        elif clade == '5. Zoosporic Fungi':
            gcas_total = len(other_fungi)
        
        prop_clade_with = gcas_with_dt/gcas_total
        st_type = 'NA'
        if 'NRPS' in st and 'PKS' in st:
            st_type = 'Multiple/Other'
        elif 'NRPS' in st:
            st_type = 'NRPS'
        elif 'PKS' in st:
            st_type = 'PKS'
        prop_pez_with = float(len(clade_subtype_gcas['1. Pezizomycotina'][st]))/len(pez_clade)
        print('\t'.join([str(x) for x in [st, prop_pez_with, clade, prop_clade_with, st_type, st]])) 

