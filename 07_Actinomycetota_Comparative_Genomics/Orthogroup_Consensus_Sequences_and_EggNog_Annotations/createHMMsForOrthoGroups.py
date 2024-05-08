import os
import sys

prot_seq_dir = os.path.abspath(sys.argv[1]) + '/'
pmsa_seq_dir = os.path.abspath(sys.argv[2]) + '/'
hmm_prof_dir = os.path.abspath(sys.argv[3]) + '/'
hmm_cons_dir = os.path.abspath(sys.argv[4]) + '/'

for f in os.listdir(prot_seq_dir):
    prot_seq_file = prot_seq_dir + f
    prot_msa_file = pmsa_seq_dir + f
    prot_hmm_file = hmm_prof_dir + f + '.hmm'
    prot_con_file = hmm_cons_dir + f
    if os.path.isfile(prot_con_file) and os.path.getsize(prot_con_file) != 0: continue
    # muscle5 super5 used for all except 4 problematic 
    align_cmd = 'muscle -super5 ' + prot_seq_file + ' -output ' + prot_msa_file + ' -threads 20 --amino' 
    #align_cmd = 'mafft --thread 20 --amino ' + prot_seq_file + ' > ' + prot_msa_file 
    hmmbuild_cmd = 'hmmbuild --amino ' + prot_hmm_file + ' ' + prot_msa_file 
    hmmemit_cmd = 'hmmemit -o ' + prot_con_file + ' -c ' + prot_hmm_file
    
    comb_cmd = align_cmd + '; ' + hmmbuild_cmd + '; ' + hmmemit_cmd
    print(comb_cmd)
