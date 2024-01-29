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
            if dom == 'PF00501.32': continue
            all_domains.add(dom)
            domain_gcas[dom].add(gca)

id_to_name = {}
with open('Pfam_Domain_Information.txt') as opdif:
    for line in opdif:
        line = line.strip()
        ls = line.split('\t')
        id_to_name[ls[0]] = ls[1]

focal_clade_ord = {'Pezizomycotina': '2', 'BGC_Enriched_Clade': '1', 'Dikarya': '2', 'Non-Dikarya Fungi': '4'}

print('Clade\tClade_Order\tDomain_ID\tDomain\tValue')
for f in os.listdir(clades_dir):
    clade_gcas = set([])
    with open(clades_dir + f) as ocf:
        for line in ocf:
            line = line.strip()
            clade_gcas.add('_'.join(line.split('_')[-2:]))
    complement_gcas = all_strains.difference(clade_gcas)
    focal_clade = f.split('.txt')[0]

    for dom in all_domains:
        dom_name = id_to_name[dom]
        focal_with = len(clade_gcas.intersection(domain_gcas[dom]))
        focal_without = len(clade_gcas.difference(domain_gcas[dom]))
        per_focal = str(focal_with/float(focal_with+focal_without))
        print(focal_clade.replace('_', ' ') + '\t' + focal_clade_ord[focal_clade] + '\t' + dom + '\t' + dom_name + '\t' + per_focal)
        if focal_clade == 'Dikarya':
            focal_with = len(complement_gcas.intersection(domain_gcas[dom]))
            focal_without = len(complement_gcas.difference(domain_gcas[dom]))
            per_focal = str(focal_with/float(focal_with+focal_without))
            print('Non-Dikarya Fungi\t4\t' + dom + '\t' + dom_name + '\t' + per_focal)
