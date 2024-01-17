import os
import sys

ufcg_results_dir = os.path.abspath('UFCG_Profile_Results/') + '/'
ufcg_tmp_dir = os.path.abspath('UFCG_Tmp_Dir/') + '/'
overview_file = 'Overview_File.txt'

with open(overview_file) as ovf:
    for i, line in enumerate(ovf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        fasta_file = ls[6]
        gca = ls[0]
        tmp_dir = ufcg_tmp_dir + gca + '/'
        os.system('mkdir ' + tmp_dir)
        print('ufcg profile -i ' + fasta_file + ' -o ' + ufcg_results_dir + gca + '/ -t 5 -w ' + tmp_dir)
