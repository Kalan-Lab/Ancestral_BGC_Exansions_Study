import os
import sys

overview_file = 'Overview_File.txt'
as_results_dir = os.path.abspath('AntiSMASH_Results/') + '/'

with open(overview_file) as ovf:
    for i, line in enumerate(ovf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0]
        fasta_file = ls[6]
        gff_file = ls[7]
        res_dir = as_results_dir + gca + '/'
        if os.path.isdir(res_dir): continue
        antismash_cmd = ['antismash', '-c', '4', '--taxon', 'fungi', '--genefinding-gff3', gff_file, '--output-dir', res_dir, fasta_file]
        print(' '.join(antismash_cmd))
