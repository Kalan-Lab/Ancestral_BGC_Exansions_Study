import os
import sys

prog = 'runFunannotateForSample.py'

all_genomes_dir = os.path.abspath(sys.argv[1]) + '/'
tmp_results_dir = os.path.abspath(sys.argv[2]) + '/'
final_results_dir = os.path.abspath(sys.argv[3]) + '/'

cpus = 5

for f in os.listdir(all_genomes_dir):
    inf = all_genomes_dir + f 
    sample = f.split('.')[0]
    samp_resdir = tmp_results_dir + sample + '/' 
    #if not os.path.isdir(final_results_dir + sample):
    print('python ' + prog + ' '  + inf + ' ' + sample + ' ' + samp_resdir + ' ' + str(cpus) + '; cp -r ' + samp_resdir + ' ' + final_results_dir + '; rm -rf ' + samp_resdir)

