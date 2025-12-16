import os
import sys
from scipy import stats

codoff_res_dir = os.path.abspath('codoff_results_v1.2.3/') + '/'
bgc_list_file = 'BGC_Details.txt'

bgc_percs = {}
for f in os.listdir(codoff_res_dir):
    bgc = f.split('.txt')[0]
    with open(codoff_res_dir + f) as ocf:
        for line in ocf:
            line = line.strip()
            ls = line.split('\t')
            if ls[0].startswith('Discordance'):
                bgc_percs[bgc] = ls[-1]

select_orders = set(['Cyanobacteriales', 'Streptosporangiales', 'Streptomycetales', 'Mycobacteriales', 'Myxococcales'])

print('BGC\tOrder\tOrder_Color\tType\tLength\tCodoff_Discordance_Percentile')
with open(bgc_list_file) as oblf:
    for i, line in enumerate(oblf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        bgc = ls[2]
        order = ls[1].split('|')[0]
        if not bgc in bgc_percs: continue
        perc = bgc_percs[bgc]
        length = ls[3]
        bgc_type = ls[-2]
        
        oc = 'Other'
        if order in select_orders:
            oc = order
        if not '|' in bgc_type:
            print(bgc + '\t' + order + '\t' + oc + '\t' + bgc_type + '\t' + str(length) + '\t' + str(perc))
        else:
            print(bgc + '\t' + order + '\t' + oc + '\tmulti-type' + '\t' + str(length) + '\t' + str(perc))
