import os
import sys
import json
from operator import itemgetter


select_file = 'Select_Genomes_for_TimeTree.txt'
ufcg_dir = 'UFCG_Profile_Results/'
prot_dir = 'Core_Protein_FASTAs/'

select = set([])
with open(select_file) as osf:
    for line in osf:
        line = line.strip()
        select.add(line)

for gca in os.listdir(ufcg_dir):
    if not gca in select: continue
    ufcg_file = None
    for f in os.listdir(ufcg_dir + gca):
        ufcg_file = ufcg_dir + gca + '/' + f
    
    assert(ufcg_file != None)
    with open(ufcg_file) as ouf:
        json_dict = json.load(ouf)
        
        for scg in json_dict['data']:
            outf = open(prot_dir + scg + '.faa', 'a+')
            all_hits = []
            for hit in json_dict['data'][scg]:
                all_hits.append([float(hit['bitscore']), hit['protein']])
            for i, hit in enumerate(sorted(all_hits, key=itemgetter(0), reverse=True)):
                if i == 0:
                    outf.write('>' + gca + '\n' + hit[1] + '\n')
            outf.close()
