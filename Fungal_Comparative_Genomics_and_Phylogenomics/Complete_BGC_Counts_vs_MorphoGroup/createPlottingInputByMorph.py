import os
import sys
from scipy import stats

as_stats_fungi = 'AntiSMASH_Stats_Fungi.txt'
as_stats_bacter = 'AntiSMASH_Stats_Bacteria.txt'
bacter_taxa_info = 'GTDB_Taxa.txt'

zyg_clade_file = 'Clades/Zygomycota.txt'
bas_clade_file = 'Clades/Basidiomycota.txt'
pez_clade_file = 'Clades/Pezizomycotina.txt'
mul_clade_file = 'Clades/Multicellular_Clade.txt'
asc_clade_file = 'Clades/Ascomycota.txt'
yea_clade_file = 'Clades/Multicellular_Clade_but_Yeast.txt'
all_fungi_file = 'Clades/All_Genomes.txt'
genus_reps_file = 'Genus_Representative_GCAs.txt'

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
yea_clade = set([])

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

with open(yea_clade_file) as oycf:
    for line in oycf:
        if line.startswith('#'): continue
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        yea_clade.add(gca)

other_fungi = all_fungi.intersection(genus_reps).difference(mul_clade)
mul_clade = mul_clade.intersection(genus_reps)
pez_clade = pez_clade.intersection(genus_reps)
asc_clade = asc_clade.intersection(genus_reps)
bas_clade = bas_clade.intersection(genus_reps)
zyg_clade = zyg_clade.intersection(genus_reps)
yea_clade = yea_clade.intersection(genus_reps)

genus_to_phylum = {}
with open(bacter_taxa_info) as obti:
    for i, line in enumerate(obti):
        line = line.strip()
        if ';g__' in line:
            genus = line.split(';g__')[1].split(';s__')[0]
            phylum = line.split(';p__')[1].split(';c__')[0]
            genus_to_phylum[genus] = phylum


print('GCA\tDomain\tPhylum_or_Clade\tComplete_BGC_Count\tBGCome_Size\tComplete_NRPS_or_PKS_BGC_Count\tColoring')
yea = []
mul = []
zoo = [] 
with open(as_stats_fungi) as oasf:
    for line in oasf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        if not gca in genus_reps: continue
        clade = 'NA'
        bgc = ls[2]
        comp_bgc = ls[3]
        bgcome_size = ls[4]
        comp_nrps_pks = ls[-1]
        if gca in pez_clade:
            color = 'Pezizomycotina'
            print(gca + '\tfungi\t1. Pezizomycotina\t' + comp_bgc + '\t' + bgcome_size + '\t' + comp_nrps_pks + '\t' + color)
        elif gca in mul_clade and not gca in pez_clade and not gca in yea_clade:
            mul.append(int(comp_bgc))
            color = 'Ascomycota'
            if gca in bas_clade:
                color = 'Basidiomycota'
            elif gca in zyg_clade:
                color = 'Zygomycota'
            print(gca + '\tfungi\t2. Multicellular (non-Pez)\t' + comp_bgc + '\t' + bgcome_size + '\t' + comp_nrps_pks + '\t' + color)
        elif gca in yea_clade:
            #yea.append(int(comp_bgc))
            color = 'Ascomycota'
            if gca in bas_clade:
                color = 'Basidiomycota'
            print(gca + '\tfungi\t3. Dikaryotic Yeast\t' + comp_bgc + '\t' + bgcome_size + '\t' + comp_nrps_pks + '\t' + color)
        elif gca in all_fungi:
            zoo.append(int(comp_bgc))
            print(gca + '\tfungi\t4. Zoosporic Fungi\t' + comp_bgc + '\t' + bgcome_size + '\t' + comp_nrps_pks + '\tZoosporic Fungi')

        if clade == 'NA': continue

print(stats.ranksums(mul, zoo, alternative='greater'))
