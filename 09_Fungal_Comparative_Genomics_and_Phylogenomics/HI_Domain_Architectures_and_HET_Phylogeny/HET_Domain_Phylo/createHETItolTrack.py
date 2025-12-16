import os
import sys
from collections import defaultdict
from ete3 import Tree
from Bio import SeqIO
import math

faa_312 = 'Proteins_with_HET_Domain.faa'
faa_ref = 'Characterized_Proteins/Characterized_HI-Associated_Proteins.faa'
het_arc = 'Domain_Architectures_of_Proteins_with_OneOrMoreHIAssociatedDoms.txt'
blast_file = 'Diamond_Blastp_to_Characterized_Proteins.txt'
pez_file = 'Clade_Info/BGC_Enriched_Pezizomycotina.txt'
aga_file = 'Clade_Info/Agaricomycetes.txt'
track_dir = 'HET_Phylo_iTol_Tracks/'

agarico = set([])
pezizo = set([])
with open(pez_file) as opf:
    for line in opf:
        line = line.strip()
        gca = ('GC' + line.split('_GC')[1].split('.')[0]).replace('GCF_', 'GCA_')
        pezizo.add(gca)

with open(aga_file) as oaf:
    for line in oaf:
        line = line.strip()
        gca = ('GC' + line.split('_GC')[1].split('.')[0]).replace('GCF_', 'GCA_')
        agarico.add(gca)

seqlen = {}
with open(faa_312) as of3:
    for rec in SeqIO.parse(of3, 'fasta'):
        seqlen[rec.id] = len(str(rec.seq))

with open(faa_ref) as ofr:
    for rec in SeqIO.parse(ofr, 'fasta'):
        seqlen[rec.id] = len(str(rec.seq))

dds = defaultdict(set)
nwp = set([])
with open(het_arc) as oha:
    for line in oha:
        line = line.strip()
        ls = line.split('\t')
        if ls[1] == 'HET|NACHT|WD40':
            nwp.add(ls[0])
        for d in ls[1].split('|'):
            dds[ls[0]].add(d)

best_hit_evalues = defaultdict(lambda: [0, 100.0])
with open(blast_file) as obf:
    for line in obf:
        line = line.strip()
        ls = line.split('\t')
        bs = float(ls[-2])
        evalue = float(ls[-3])
        if bs > best_hit_evalues[ls[0]][0]:
            best_hit_evalues[ls[0]] = [bs, evalue]
        elif bs == best_hit_evalues[ls[0]][0]:
            best_hit_evalues[ls[0]] = [bs, min([best_hit_evalues[ls[0]][1], evalue])]

outf1 = open(track_dir + 'Distinct_Domains.iTol.txt', 'w')
outf2 = open(track_dir + 'Protein_Lengths.iTol.txt', 'w')
outf3 = open(track_dir + 'Pod_HET-E_Structure.iTol.txt', 'w')
outf4 = open(track_dir + 'DBlastp_Evalues.iTol.txt', 'w')
outf5 = open(track_dir + 'Characterized_Piecharts.iTol.txt', 'w')
outf6 = open(track_dir + 'Taxonomy.iTol.txt', 'w')

outf1.write('DATASET_SIMPLEBAR\n')
outf1.write('SEPARATOR TAB\n')
outf1.write('COLOR\t#000000\n')
outf1.write('DATASET_LABEL\tDistinct Pfam Domains\n')
outf1.write('DATA\n')

outf2.write('DATASET_SIMPLEBAR\n')
outf2.write('SEPARATOR TAB\n')
outf2.write('COLOR\t#000000\n')
outf2.write('DATASET_LABEL\tProtein Length\n')
outf2.write('DATA\n')

outf3.write('DATASET_COLORSTRIP\n')
outf3.write('SEPARATOR TAB\n')
outf3.write('COLOR\t#000000\n')
outf3.write('DATASET_LABEL\tHET-NACHT-WD40\n')
outf3.write('DATA\n')

outf4.write('DATASET_HEATMAP\n')
outf4.write('SEPARATOR TAB\n')
outf4.write('COLOR\t#000000\n')
outf4.write('DATASET_LABEL\tDIAMOND Blastp E-value\n')
outf4.write('FIELD_LABELS\tE-value\n')
outf4.write('DATA\n')

outf5.write('DATASET_PIECHART\n')
outf5.write('SEPARATOR TAB\n')
outf5.write('COLOR\t#000000\n')
outf5.write('DATASET_LABEL\tRef\n')
outf5.write('FIELD_LABELS\tRef\n')
outf5.write('FIELD_COLORS\t#eb4034\n')
outf5.write('DATA\n')

outf6.write('DATASET_COLORSTRIP\n')
outf6.write('SEPARATOR TAB\n')
outf6.write('COLOR\t#000000\n')
outf6.write('DATASET_LABEL\tTaxonomic_Clade\n')
outf6.write('DATA\n')

t = Tree('MAFFT_add_to_MUSCLE_align.tre')
for n in t.traverse('postorder'):
    if not n.is_leaf(): continue
    dn = n.name
    pn = dn
    is_ref = "1"
    gca = None
    if 'CDS' in dn:
        pn = '|'.join(dn.split('|')[:-2])
        is_ref = "0"
        gca = pn.split('|')[0]

    distinct_domains = len(dds[pn])
    protein_length = seqlen[pn]
    
    da_type = 'other'
    da_type_color = '#FFFFFF'
    if pn in nwp:
        da_type = 'het-nacht-wd40'
        da_type_color = '#000000'

    outf1.write(dn + '\t' + str(distinct_domains) + '\n')
    outf2.write(dn + '\t' + str(protein_length) + '\n')
    outf3.write(dn + '\t' + da_type_color + '\t' + da_type + '\n')
    outf4.write(dn + '\t' + str(0-math.log(best_hit_evalues[dn][1], 10)) + '\n')
    if is_ref == '1':
        outf5.write(dn + '\t1\t1\t' + is_ref + '\n')
    if gca in pezizo:
        outf6.write(dn + '\t#A8CCF4\tbgc_rich_pezizo\n')
    elif gca in agarico:
        outf6.write(dn + '\t#EB6964\tagaricomycetes\n')
    else:
        print(gca)

outf1.close()
outf2.close()
outf3.close()
outf4.close()
outf5.close()
outf6.close()
