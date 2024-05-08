import os
import sys
from collections import defaultdict
import statistics

resfams_dir = 'ResFam_Search_Results/' # directory not included with HMMER3 results for each of the ~33k actino proteome files from searching with resfams core HMM set 
gca_to_taxa_file = '../GCA_to_Taxonomy.txt'

all_samples = set([])
all_rfs = set([])
sample_rf_lts = defaultdict(lambda: defaultdict(set))
for f in os.listdir(resfams_dir):
    s = f.split('.')[0]
    all_samples.add(s)
    with open(resfams_dir + f) as ordf:
        for line in ordf:
            line = line.strip()
            if line.startswith('#'): continue
            ls = line.split()
            lt = ls[0]
            rf = ls[3]
            sample_rf_lts[s][rf].add(lt)
            all_rfs.add(rf)

gca_to_taxa = {}
with open(gca_to_taxa_file) as ogtf:
    for line in ogtf:
        line = line.strip()
        ls = line.split('\t')
        gca_to_taxa[ls[0].split('.')[0]] = ls[1].split(';f__')[1].split(';g__')[0]

family_counts = defaultdict(lambda: defaultdict(list))
for s in all_samples:
    f = gca_to_taxa[s]
    for rf in all_rfs:
        srfc = len(sample_rf_lts[s][rf])
        family_counts[f][rf].append(srfc)

print('DATASET_HEATMAP')
print('SEPARATOR TAB')
print('DATASET_LABEL\tResFams')
print('COLOR\t#000000')
print('FIELD_LABELS\t' + '\t'.join(sorted(list(all_rfs))))
print('DATA')
for f in family_counts:
    printlist = [f]
    for rf in sorted(list(all_rfs)):
        printlist.append(str(statistics.mean(family_counts[f][rf])))
    print('\t'.join(printlist))
