import os
import sys

genomad_dir = '/workspace/local/rauf/multicellularity_and_bgcs/genomad/geNomad_Results/'

print('gca\tphage_sum\tplasmid_sum\ttotal_sum')
for s in os.listdir(genomad_dir):
    samp_dir = genomad_dir + s  + '/'
    if os.path.isdir(samp_dir):
        gca = s.split('.')[0]
        plas_file = samp_dir + s + '_summary/' + s + '_plasmid_summary.tsv'
        phage_file = samp_dir + s + '_summary/' + s + '_virus_summary.tsv'
        assert(os.path.isfile(plas_file) and os.path.isfile(phage_file))
        phage_sum = 0
        plasmid_sum = 0
        plasmid_scaffs = set([])
        total_sum = 0
        with open(plas_file) as opf:
            for i, line in enumerate(opf):
                line = line.strip()
                ls = line.split('\t')
                if i == 0: continue
                plasmid_scaffs.add(ls[0])
                plasmid_sum += int(ls[1])
                total_sum += int(ls[1])

        with open(phage_file) as orf:
            for i, line in enumerate(orf):
                line = line.strip()
                ls = line.split('\t')
                if i == 0: continue
                scaff = ls[0].split('|')[0]
                if not scaff in plasmid_scaffs:
                    total_sum += int(ls[1])
                phage_sum += int(ls[1])
            print(gca + '\t' + str(phage_sum) + '\t'  + str(plasmid_sum) + '\t' + str(total_sum))
