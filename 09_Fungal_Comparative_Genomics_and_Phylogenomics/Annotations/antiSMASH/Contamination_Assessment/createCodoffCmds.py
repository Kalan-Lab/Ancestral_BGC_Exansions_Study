import os
import sys
from collections import defaultdict

as_res_dir = os.path.abspath('../All_antiSMASH_Results/') + '/'
res_dir = os.path.abspath('Codoff_Results/') + '/'
for s in os.listdir(as_res_dir): 
    samp_dir = as_res_dir + s + '/'
    samp_res_dir = res_dir + s + '/'
    os.system('mkdir ' + samp_res_dir)
    if os.path.isdir(samp_dir):
        full_genome_gbk = None
        for f in os.listdir(samp_dir):
            if f.endswith('.gbk') and not 'region' in f:
                full_genome_gbk = samp_dir + f
        for f in os.listdir(samp_dir):
            if f.endswith('.gbk') and 'region' in f:
                region_gbk = samp_dir + f 
                res_file = samp_res_dir + f
                print('codoff -g ' + full_genome_gbk + ' -f ' + region_gbk + ' -o ' + res_file) 
