import os
import sys
from Bio import SeqIO

msa_dir = os.path.abspath('zol_dom_Results/Ortholog_Group_Processing/OG_Protein_Alignments_Trimmed/') + '/'
res_dir = os.path.abspath('sc_cons_msas/') + '/'
select_file = 'sc_domains_90per_conserved.txt'

select = set([])
with open(select_file) as osf:
    for line in osf:
        line = line.strip()
        select.add(line.split('\t')[0][1:])

for f in os.listdir(msa_dir):
    dog = f.split('.')[0]
    if dog not in select: continue
    outf = open(res_dir + f, 'w')
    with open(msa_dir + f) as omf:
        for rec in SeqIO.parse(omf, 'fasta'):
            outf.write('>' + rec.id.split('|')[0].split('.')[0] + '\n' + str(rec.seq) + '\n')
    outf.close() 
