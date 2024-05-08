import os
import sys
from Bio import SeqIO

seq_algn_dir = os.path.abspath(sys.argv[1]) + '/'
seq_filt_algn_dir = os.path.abspath(sys.argv[2]) + '/'

for f in os.listdir(seq_algn_dir):
    os.system('trimal -in ' + seq_algn_dir + f + ' -out ' + seq_filt_algn_dir + f + ' -strict -keepseqs')  
