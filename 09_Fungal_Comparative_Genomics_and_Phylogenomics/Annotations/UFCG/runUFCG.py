import os
import sys

ufcg_results_dir = os.path.abspath('UFCG_Profile_Results/') + '/'
ufcg_tmp_dir = os.path.abspath('UFCG_Tmp_Dir/') + '/'
local_input_dir = os.path.abspath('Local_Inputs/') + '/'
genome_listing_file = 'Overview_File.txt'

with open(genome_listing_file) as oglf:
    for i, line in enumerate(oglf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        fasta_file = ls[6]
        gca = ls[0]
        tmp_dir = ufcg_tmp_dir + gca + '/'
        res_dir = ufcg_results_dir + gca + '/'
        if not os.path.isdir(res_dir):
            assert('/local/' in fasta_file)
            os.system('rm -rf ' + tmp_dir)
            os.system('mkdir ' + tmp_dir)
            fasta_file = local_input_dir + fasta_file.split('/')[-1]
            print('ufcg profile -i ' + fasta_file + ' -o ' + ufcg_results_dir + gca + '/ -t 4 -w ' + tmp_dir)
