import os
import sys

tiny_data_file = 'Tiny_AAI_for_Plotting.txt'
gbk_dir = 'fai_Results/Final_Results/Homologous_Gene_Cluster_GenBanks//'
res_dir = 'zol_input/'

select = set([])
gca_to_clade = {}
with open(tiny_data_file) as otdf:
    for i, line in enumerate(otdf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca = ls[1].split('.')[0]
        aai = float(ls[2])
        sp = float(ls[3])
        if sp >= 0.25 and aai >= 60.0:
            select.add(gca)
            gca_to_clade[gca] = ls[0]

for f in os.listdir(gbk_dir):
    gca = f.split('.')[0]
    gbk_file = gbk_dir + f
    if gca in select:
        if gca_to_clade[gca] == 'Actinomycetia - Clade-2':
            print(f)
        os.system('cp ' + gbk_file + ' ' + res_dir)
