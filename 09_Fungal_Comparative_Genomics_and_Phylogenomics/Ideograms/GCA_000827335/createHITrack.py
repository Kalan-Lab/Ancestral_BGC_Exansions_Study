import os
import sys
from Bio import SeqIO

prot_locs = {}
with open('GCA_000827335.faa') as opf:
    for rec in SeqIO.parse(opf, 'fasta'):
        desc = rec.description
        prot = rec.id
        chrom = desc.split()[1]
        midpoint = ((int(desc.split()[2])-1) + (int(desc.split()[3])-1))/2.0
        prot_locs[prot] = [chrom, midpoint]

print('\t'.join(['Protein', 'Chrom', 'Position', 'Type']))

dom_architect = set(['HET', 'Het-C', 'ATP-cone|Ribonuc_red_lgN|Ribonuc_red_lgC', 'HMG_box', 'MAT1-1-2', 'Patatin|NB-ARC', 'HeLo|HET-s_218-289', 'HET|NACHT|WD40', 'SET|Rubis-subs-bind|zf-MYND', 'Peptidase_S8', 'GLTP', 'RVT_2', 'SERF-like_N', 'CYSTM', 'Ccdc124', 'Fructosamin_kin', 'Patatin|TPR_10'])

print('\t'.join(['Protein', 'Chrom', 'Position', 'Type']))
with gzip.open("Domain_Architectures_of_Proteins_with_OneOrMoreHIAssociatedDoms.txt.gz", 'rt') as odf:
    for line in odf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('|')[0]
        if gca == 'GCA_000827335':
            prot = ls[0]
            dom = ls[1]
            if dom in dom_architect:
                dtype = 'HET'
                if dom != 'HET':
                    dtype = 'Other'
                chrom, pos = prot_locs[prot]
                print('\t'.join([prot, chrom, str(pos), dtype]))

"""
with open('HIT_Hits.Domain.Cut_TC.txt') as ohf:
    for line in ohf:
        if line.startswith('#'): continue
        line = line.strip()
        ls = line.split()
        hit_type = ls[3]
        prot = ls[0]
        chrom, pos = prot_locs[prot]
        print('\t'.join([prot, chrom, str(pos), hit_type]))
"""
