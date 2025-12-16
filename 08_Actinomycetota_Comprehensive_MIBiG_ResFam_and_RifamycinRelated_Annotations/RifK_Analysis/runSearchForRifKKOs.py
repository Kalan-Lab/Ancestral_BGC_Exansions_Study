import os
import sys

prot_dir = os.path.abspath(sys.argv[1]) + '/'
results_dir = os.path.abspath(sys.argv[2]) + '/'
tmp_dir = os.path.abspath(sys.argv[3]) + '/'

for f in os.listdir(prot_dir):
    inf = prot_dir + f
    outf = results_dir + f + '.txt'
    if os.path.isfile(outf): continue
    tmp_dir_sample = tmp_dir + f + '/'
    os.system('mkdir ' + tmp_dir_sample)
    cmd = 'exec_annotation -p KOFam_RifK_HMM/ -k ko_list --cpu 1 --tmp-dir ' + tmp_dir_sample + ' -o ' + outf + ' ' + inf
    print(cmd)
