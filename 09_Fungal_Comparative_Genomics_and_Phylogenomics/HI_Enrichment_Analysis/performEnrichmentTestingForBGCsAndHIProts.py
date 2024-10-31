import os
import sys
from Bio import SeqIO
from scipy import stats

mc_gcas = set([])
with open('Agaricomycetes.txt') as oaf:
    for line in oaf:
        gca = '_'.join(line.strip().split('.')[0].split('_')[-2:])
        mc_gcas.add(gca)

with open('BGC_Enriched_Pezizomycotina.txt') as oaf:
    for line in oaf:
        gca = '_'.join(line.strip().split('.')[0].split('_')[-2:])
        mc_gcas.add(gca)

bgc_proteins = set([])
with open('BGC_Proteins.txt') as of:
    for line in of:
        line = line.strip()
        gca = line.split('|')[0]
        if gca in mc_gcas:
            bgc_proteins.add(line)

other_proteins = set([])
with open('../All_Proteins.faa') as of:
    for rec in SeqIO.parse(of, 'fasta'):
        prot = rec.id
        gca = prot.split('|')[0]
        if gca in mc_gcas and not prot in bgc_proteins:
            other_proteins.add(prot)

bgc_hi = 0
other_hi = 0
with open('HIT_Hits.Domain.Cut_TC.txt') as of:
    for line in of:
        if line.startswith('#'): continue
        line = line.strip()
        prot = line.split()[0]
        if prot in bgc_proteins:
            bgc_hi += 1
        elif prot in other_proteins:
            other_hi += 1

bgc_other = len(bgc_proteins) - bgc_hi
other_other = len(other_proteins) - other_hi

print('BGC HI\t' + str(bgc_hi))
print('BGC Other\t' + str(bgc_other))
print('Other HI\t' + str(other_hi))
print('Other Other\t' + str(other_other))

print(stats.fisher_exact([[bgc_hi, bgc_other], [other_hi, other_other]], alternative='greater'))
