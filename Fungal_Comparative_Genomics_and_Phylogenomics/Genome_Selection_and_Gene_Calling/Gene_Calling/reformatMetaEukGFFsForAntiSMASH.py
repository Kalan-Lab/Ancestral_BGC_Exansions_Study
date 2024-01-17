import os
import sys
from collections import defaultdict

metaeuk_results_dir = os.path.abspath(sys.argv[1]) + '/' #'NonAscomycota_MetaEuk_Results/'
process_gff_dir = os.path.abspath(sys.argv[2]) + '/' #'AntiSMASH_GFF_Inputs/'

for f in os.listdir(metaeuk_results_dir):
    if not f.endswith('.gff'): continue
    metaeuk_results_file = metaeuk_results_dir + f
    gff_processed_file = process_gff_dir + f 
    gpf_outf = open(gff_processed_file, 'w')
    coord_data = defaultdict(lambda: defaultdict(dict))
    with open(metaeuk_results_file) as omrf:
        for line in omrf:
            line = line.strip()
            if line.startswith('#'): 
                gpf_outf.write(line + '\n')
            else:
                ls = line.split('\t')
                print(ls)
                feature_type = ls[2]
                if not feature_type == 'CDS': continue
                scaffold = ls[0]
                start = ls[3]
                end = ls[4]
                score = float(ls[5])
                coord = tuple([scaffold, start, end])
                coord_data[coord][score] = line

    for coord in coord_data:
        for i, score in enumerate(sorted(coord_data[coord], reverse=True)):
            if i > 0: continue
            gpf_outf.write(coord_data[coord][score] + '\n')
    gpf_outf.close()
