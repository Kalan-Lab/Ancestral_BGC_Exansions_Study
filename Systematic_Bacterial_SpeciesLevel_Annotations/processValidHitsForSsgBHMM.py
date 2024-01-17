import os
import sys

hmm_results_dir = 'Search_Results/'
score_threshold = 99.9

for f in os.listdir(hmm_results_dir):
    gca = f.split('.')[0]
    with open(hmm_results_dir + f) as of:
        for line in of:
            line = line.strip()
            if line.startswith('#'): continue
            ls = line.split()
            score = float(ls[5])
            if score >= score_threshold:
                print(gca)
