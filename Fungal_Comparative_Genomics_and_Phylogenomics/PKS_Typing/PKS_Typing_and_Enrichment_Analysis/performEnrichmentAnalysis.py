import os
import sys
from collections import defaultdict
from scipy import stats

all_strains_file = 'All_Genomes.txt'
clades_dir = 'Clades/'
domain_dir = 'Results_Domain_Presence/'

all_strains = set([])
with open(all_strains_file) as oasf:
    for line in oasf:
        line = line.strip()
        all_strains.add(line)

all_domains = set([])
domain_gcas = defaultdict(set)
for gca in os.listdir(domain_dir):
    dom_file = domain_dir + gca + '/domains_found.txt'
    if not os.path.isfile(dom_file): continue
    with open(dom_file) as odf:
        for line in odf:
            dom = line.strip() 
            all_domains.add(dom)
            domain_gcas[dom].add(gca)

for f in os.listdir(clades_dir):
    clade_gcas = set([])
    with open(clades_dir + f) as ocf:
        for line in ocf:
            line = line.strip()
            clade_gcas.add('_'.join(line.split('_')[-2:]))
    complement_gcas = all_strains.difference(clade_gcas)
    focal_clade = f.split('.txt')[0]
    for dom in all_domains:
        focal_with = len(clade_gcas.intersection(domain_gcas[dom]))
        focal_without = len(clade_gcas.difference(domain_gcas[dom]))
        other_with = len(complement_gcas.intersection(domain_gcas[dom]))
        other_without = len(complement_gcas.difference(domain_gcas[dom]))
        _, pvalue = stats.fisher_exact([[focal_with, focal_without], [other_with, other_without]])
        per_focal = str(focal_with/float(focal_with+focal_without))
        per_other = str(other_with/float(other_with+other_without))
        print(focal_clade + '\t' + dom + '\t' + str(pvalue) + '\t' + per_focal + '\t' + per_other)
