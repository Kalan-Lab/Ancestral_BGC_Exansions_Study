import os
import sys

rae_query_file = 'RAE_Palindrome.fna'
actino_genome_dir = os.path.abspath('All_Actino_Genomes') + '/'
blast_res_dir = os.path.abspath('BLASTn_Results') + '/'

for f in os.listdir(actino_genome_dir):
    print('blastn -query ' + rae_query_file + ' -subject ' + actino_genome_dir + f + ' -out ' + blast_res_dir + f + '.txt' + ' -outfmt "6 std qcovhsp" -ungapped -task blastn-short')
