import os
import sys
from Bio import SeqIO

prog = 'determinePKSType.py'
res_dir = 'Results/'
hmm_db = os.path.abspath('Fungal_PKS_Typing.hmms')

antismash_dir = os.path.abspath('../AntiSMASH_Results/') + '/'
for rec in SeqIO.parse('../KS_Proteins.faa', 'fasta'):
    s = rec.id.split('|')[0]
    reg = rec.id.split('|')[1]
    bgc_gbk = antismash_dir + s + '/' + reg + '.gbk'
    #os.system('mkdir ' + res_dir + s + '_bgc_' + reg + '/')
    print(bgc_gbk)
    #os.system('python ' + prog + ' ' + bgc_gbk + ' ' + hmm_db + ' ' + res_dir + s + '_bgc_' + reg + '/')
