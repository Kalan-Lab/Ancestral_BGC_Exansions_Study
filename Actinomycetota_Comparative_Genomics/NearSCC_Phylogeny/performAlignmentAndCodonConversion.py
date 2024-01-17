import os
import sys
from Bio import SeqIO
from collections import defaultdict

core_ogs_file = 'NearSingleCopyCore_OGs.txt'
og_dir = 'Orthogroup_Sequences/'
og_msa_dir = 'Orthogroup_Protein_Alignments/'

core_ogs = set([])
with open(core_ogs_file) as ocof:
    for line in ocof:
        line = line.strip()
        ls = line.split('\t')
        core_ogs.add(ls[0])

print(core_ogs)
og_to_prot = defaultdict(dict)
samples = set([])
for f in os.listdir(og_dir):
    og = f.split('.fa')[0]
    if not og in core_ogs: continue
    with open(og_dir + f) as ogf:
        for rec in SeqIO.parse(ogf, 'fasta'):
            samp = rec.id.split('|')[0]
            og_to_prot[og][samp] = rec.id
            samples.add(samp)
    pro_og_msa = og_msa_dir + f
    muscle_cmd = 'muscle -align %s -output %s -threads 30' % (og_dir + f, pro_og_msa)
    os.system(muscle_cmd)
