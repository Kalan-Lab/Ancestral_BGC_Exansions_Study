import os
import sys
from Bio import SeqIO

prot_dir = 'Core_Protein_FASTAs/'
prot_msa_dir = 'Core_Protein_MSAs/'
prot_trim_dir = 'Core_Protein_MSAs_Trimmed/'

for f in os.listdir(prot_dir):
    prot_file = prot_dir + f

    recs = 0
    with open(prot_file) as opf:
        for rec in SeqIO.parse(opf, 'fasta'):
            recs += 1 
    if recs < 285: continue

    msa_file = prot_msa_dir + f 
    trim_file = prot_trim_dir + f

    muscle_cmd = ['muscle', '-align', prot_file, '-output', msa_file, '-threads', '40']
    trimal_cmd = ['trimal', '-in', msa_file, '-out', trim_file, '-keepseqs', '-strict'] # was originally using 90% un-ambiguous cutoff
    os.system(' '.join(muscle_cmd) + '; ' + ' '.join(trimal_cmd))
