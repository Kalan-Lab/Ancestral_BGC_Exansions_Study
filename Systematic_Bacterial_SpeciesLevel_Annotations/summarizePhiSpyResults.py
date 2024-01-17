import os
import sys

phispy_dir = 'PhiSpy_Results_21k/'

for s in os.listdir(phispy_dir):
    samp_dir = phispy_dir + s  + '/'
    if os.path.isdir(samp_dir):
        gca = s.split('.')[0]
        res_file = samp_dir + 'prophage_coordinates.tsv'
        if os.path.isfile(res_file):
            sum_phage = 0
            phage_count = 0
            with open(res_file) as orf:
                for line in orf:
                    line = line.strip()
                    ls = line.split('\t')
                    sum_phage += abs(int(ls[3])-int(ls[2]))
                    phage_count += 1

            print(gca + '\t' + str(phage_count) + '\t' + str(sum_phage))
        else:
            print(gca + '\t0\t0')
