import os
import sys
from Bio import SeqIO

prog = 'determineGenomeWidePresenceOfRelevantDomains.py'
res_dir = 'Results_GenomeWide_Domain_Presence/'
hmm_db = os.path.abspath('Fungal_PKS_Typing.hmms')

antismash_dir = os.path.abspath('../AntiSMASH_Results/') + '/'
for s in os.listdir(antismash_dir):
    samp_dir = antismash_dir + s + '/'
    for f in os.listdir(samp_dir):
        if f.endswith('.gbk') and not '.region' in f:
            gbk_file = samp_dir + f
            os.system('mkdir ' + res_dir + s + '/')
            os.system('python ' + prog + ' ' + gbk_file + ' ' + hmm_db + ' ' + res_dir + s + '/')
