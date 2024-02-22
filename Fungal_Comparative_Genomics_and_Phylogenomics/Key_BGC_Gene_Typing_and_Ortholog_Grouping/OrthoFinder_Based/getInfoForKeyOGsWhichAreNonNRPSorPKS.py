import os
import sys
from collections import defaultdict
from operator import itemgetter

pks_or_nrps_ogs = set(['OG0000000', 'OG0000002', 'OG0000004', 'OG0000006'])

bas_clade_file = 'Taxonomic_Clades/Basidiomycota.txt'
zyg_clade_file = 'Taxonomic_Clades/Zygomycota.txt'
pez_clade_file = 'Taxonomic_Clades/Pezizomycotina.txt'
mul_clade_file = 'Taxonomic_Clades/Multicellular_Clade.txt'
asc_clade_file = 'Taxonomic_Clades/Ascomycota.txt'
all_fungi_file = 'Taxonomic_Clades/All_Genomes.txt'
genus_reps_file = 'Genus_Representative_GCAs.txt'

orthogroup_file = 'OrthoFinder_Results/Results_Feb12/Orthogroups/Orthogroups.tsv'
domain_annot_file = 'Best_DomainMatches_for_Common_OGs.txt'
og_annot_file = 'Select_OG_Annots.txt'

genus_reps = set([])
with open(genus_reps_file) as ogrf:
    for line in ogrf:
        line = line.strip()
        genus_reps.add(line)

mul_clade = set([])
pez_clade = set([])
all_fungi = set([])
bas_clade = set([])
zyg_clade = set([])
asc_clade = set([])

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

og_annot = {}
with open(og_annot_file) as oaf:
    for line in oaf:
        line = line.strip()
        ls = line.split('\t')
        og_annot[ls[0]] = ls[1]


print('OG\tOG_Order\tClade\tClade_Proportion_With\tBGC_Class\tOG_Name')
for line in open(orthogroup_file):
    if line.startswith('Orthogroup'): continue
    line = line.strip('\n')
    ls = line.split('\t')
    og = ls[0]
    gc = 0
    clade_count = defaultdict(int)
    type_count = defaultdict(int)
    type_class_count = defaultdict(int)
    for gs in ls[1:]:
        if gs.strip() != '':
            gc += 1
            gca = gs.split('|')[0]
            if not gca in genus_reps: continue
            for g in gs.split(', '):
                for t in g.split('|')[3:]:
                    tc = 'Multiple/Other'
                    if 'NRP' in t:
                        tc = 'NRPS'
                    elif 'PKS' in t:
                        tc = 'PKS'
                    elif 'terpene' == t:
                        tc = 'Terpene'
                    type_count[t] += 1
                    type_class_count[tc] += 1
            if gca in pez_clade:
                clade_count['pez'] += 1.0/len(pez_clade)
            elif gca in asc_clade:
                clade_count['asc'] += 1.0/len(asc_clade)
            elif gca in bas_clade:
                clade_count['bas'] += 1.0/len(bas_clade)
            elif gca in zyg_clade:
                clade_count['zyg'] += 1.0/len(zyg_clade)
            elif gca in other_fungi:
                clade_count['zoo'] += 1.0/len(other_fungi)

    type_string = []
    for t in sorted(type_count.items(), key=itemgetter(1), reverse=True):
        type_string.append(t[0] + '=' + str(t[1]))

    top_tc = "Multiple/Other"
    tot = sum(type_class_count.values())
    for tc in sorted(type_class_count):
        tc_prop = type_class_count[tc]/float(tot)
        if tc_prop >= 0.7:
            top_tc = tc

    if og in pks_or_nrps_ogs: continue
    count_string = []
    for c in sorted(clade_count):
        count_string.append(c + '=' + str(clade_count[c]))
    if clade_count['pez'] >= 0.50:
        print(og + '\t' + str(clade_count['pez']) + '\t1. Pezizomycotina\t' + str(clade_count['pez']) + '\t' + top_tc + '\t' + og_annot[og])
        print(og + '\t' + str(clade_count['pez']) + '\t2. Other Ascomycota\t' + str(clade_count['asc']) + '\t' + top_tc + '\t' + og_annot[og])
        print(og + '\t' + str(clade_count['pez']) + '\t3. Basidiomycota\t' + str(clade_count['bas']) + '\t' + top_tc + '\t' + og_annot[og])
        print(og + '\t' + str(clade_count['pez']) + '\t4. Zygomycota\t' + str(clade_count['zyg']) + '\t' + top_tc + '\t' + og_annot[og])
        print(og + '\t' + str(clade_count['pez']) + '\t5. Zoosporic Fungi\t' + str(clade_count['zoo']) + '\t' + top_tc + '\t' + og_annot[og])
