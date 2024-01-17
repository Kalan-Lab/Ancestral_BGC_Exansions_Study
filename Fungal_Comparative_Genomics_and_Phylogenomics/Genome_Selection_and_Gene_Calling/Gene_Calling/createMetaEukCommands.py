import os
import sys

genome_dir = 'Genomes/' # Directory of input genomes in FASTA format
me_db_file = 'MetaEuk/database/MMETSP_uniclust50_MERC'
tmp_results_dir = 'MetaEuk_Tmp_Results/'
fin_results_dir = 'MetaEuk_Results/'

genome_paths = {}
for f in os.listdir(genome_dir):
    genome_file = genome_dir + f
    gca = '_'.join(f.split('_')[:2])
    res_file = tmp_results_dir + gca + '_RESULT_'
    check_res_file = fin_results_dir + gca + '_RESULT_.gff'
    tmp_dir = tmp_results_dir + gca + '_TMPDIR/'
    os.system('mkdir ' + tmp_dir)
    if not os.path.isfile(check_res_file):
        cmd = ['metaeuk', 'easy-predict', '--threads', '20', genome_file, me_db_file, res_file, tmp_dir, ';', 'rm', '-rf', tmp_dir, ';', 'cp', res_file + '*', fin_results_dir]
        print(' '.join(cmd))
