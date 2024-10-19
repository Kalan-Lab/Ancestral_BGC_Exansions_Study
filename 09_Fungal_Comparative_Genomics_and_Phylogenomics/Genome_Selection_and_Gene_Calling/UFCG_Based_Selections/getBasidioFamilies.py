import os
import sys

ufcg_dir = 'Species_UFCG/'
for f in os.listdir(ufcg_dir):
    gca = '_'.join(f.split('_')[:2])
    with open(ufcg_dir + f) as ocf:
        for line in ocf:
            line = line.strip()
            if line.startswith('"taxonomy": ') and 'Basidiomycota' in line and not 'fam_Incertae_sedis' in line and not 'unidentified' in line:
                print(line.split(';')[4])

