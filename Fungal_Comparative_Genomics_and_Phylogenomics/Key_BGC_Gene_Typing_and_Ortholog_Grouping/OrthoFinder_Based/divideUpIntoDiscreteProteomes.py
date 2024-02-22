import os
import sys
from Bio import SeqIO

dom_faa = 'All_Key_Proteins.faa'
prot_dir = os.path.abspath('BGComes') + '/'

with open(dom_faa) as odf:
    for rec in SeqIO.parse(odf, 'fasta'):
        gca = rec.id.split('|')[0]
        prot_file = prot_dir + gca + '.faa'
        outf = open(prot_file, 'a+')
        outf.write('>' + rec.id + '\n' + str(rec.seq) + '\n')
        outf.close()
