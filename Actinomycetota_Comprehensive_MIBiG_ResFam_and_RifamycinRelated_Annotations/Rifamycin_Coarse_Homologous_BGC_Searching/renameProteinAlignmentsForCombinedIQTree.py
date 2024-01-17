import os
import sys
from Bio import SeqIO

select_ogs_file = 'select_ogs.txt'
gca_to_tax_file = '../GCA_to_Taxonomy.txt'
tiny_plot_file = 'Tiny_AAI_for_Plotting.txt'
prot_dir = os.path.abspath('zol_results/Ortholog_Group_Processing/OG_Protein_Alignments_Trimmed/') + '/'
rename_dir = os.path.abspath('Combined_IQ_Tree_Input/') + '/'

select = set([])
with open(select_ogs_file) as osof:
    for line in osof:
        line = line.strip()
        select.add(line)

gca_to_genus = {}
with open(gca_to_tax_file) as ogttf:
    for line in ogttf:
        line = line.strip()
        ls = line.split('\t')
        gca_to_genus[ls[0]] = ls[1].split(';g__')[1].split(';s__')[0]

gca_to_clade = {}
with open(tiny_plot_file) as otpf:
    for line in otpf:
        line = line.strip()
        ls = line.split('\t')
        if len(ls[0].split(' - ')) >= 2:
            gca_to_clade[ls[1].split('.')[0]] = ls[0].split(' - ')[1]

for f in os.listdir(prot_dir):
    og = f.split('.msa.faa')[0]
    if not og in select: continue
    outf = open(rename_dir + f, 'w')
    with open(prot_dir + f) as opf:
        for rec in SeqIO.parse(opf, 'fasta'):
            gca = rec.id.split('.')[0]
            genus = gca_to_genus[gca]
            clade = gca_to_clade[gca]
            outf.write('>' + gca + '|' + genus + '|' + clade + '\n' + str(rec.seq) + '\n')
    outf.close()
