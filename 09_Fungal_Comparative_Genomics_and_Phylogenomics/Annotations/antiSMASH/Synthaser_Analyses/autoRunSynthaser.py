import os
import sys

input_dir =  os.path.abspath('synthaser_inputs/') + '/'
result_dir = os.path.abspath('synthaser_results/') + '/'

for f in os.listdir(input_dir):
    sample = f.split('.faa')[0]
    inf = input_dir + f
    outf = result_dir + sample + '.txt'
    json = result_dir + sample + '.json'

    os.system('synthaser search -json ' + json + ' -o ' + outf + ' -qf ' + inf)
