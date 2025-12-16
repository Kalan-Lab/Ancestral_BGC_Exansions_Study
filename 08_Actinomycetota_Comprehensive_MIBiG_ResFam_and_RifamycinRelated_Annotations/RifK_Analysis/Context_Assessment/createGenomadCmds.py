import os
import sys

genomad_db = '/workspace/local/rauf/multicellularity_and_bgcs/genomad/genomad_db/'
results_dir = 'geNomad_Results/'


for f in os.listdir(genome_dir):
    if not f.endswith('.fna'): continue
    input_genome = genome_dir + f
    genomad_results = results_dir + '.'.join(f.split('.')[:-1]) + '/'

    genomad_cmd = ['genomad', 'end-to-end', '--cleanup', '--threads', '4', '--splits', '8', input_genome, genomad_results, genomad_db]
    print(' '.join(genomad_cmd))

