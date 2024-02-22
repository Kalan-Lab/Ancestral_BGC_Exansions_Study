import os
import sys
from collections import defaultdict
from Bio import SeqIO
import json 

pez_clade_file = 'Clades/Pezizomycotina.txt'
all_fungi_file = 'Clades/All_Genomes.txt'
all_prots_file = '../All_PKS_Full_Proteins.faa'

pez_clade = set([])
all_fungi = set([])

with open(all_fungi_file) as oaff:
    for line in oaff:
        line = line.strip()
        all_fungi.add(line)

with open(pez_clade_file) as obef:
    for line in obef:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        pez_clade.add(gca)

other_fungi = all_fungi.difference(pez_clade)

id_to_gca = {}
with open(all_prots_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        id_to_gca[rec.id] = rec.description.split('|')[1]

synth_files = ['synthaser_calls_on_pks_p1.json', 'synthaser_calls_on_pks_p2.json', 'synthaser_calls_on_pks_p3.json', 'synthaser_calls_on_pks_p4.json']


clade_domain_gcas = defaultdict(lambda: defaultdict(set))
all_dts = set([])
for sf in synth_files:
    with open(sf) as osf:
        d = json.load(osf)
        for pks in d:
            pks_id = pks['header']
            pks_gca = id_to_gca[pks_id]
            pks_clade = None
            pks_seq = pks['sequence']
            if pks_gca in pez_clade:
                pks_clade = 'Pezizomycotina'
            elif pks_gca in other_fungi:
                pks_clade = 'Other'
            if pks_clade != None and len(pks['classification']) > 0 and pks['classification'][0] == 'PKS' and 'Non-reducing' in pks['classification']:
                for dom in pks['domains']:
                    dt = dom['type']
                    if dt == 'PT':
                        start = dom['start']
                        end = dom['end']
                        dom_seq = pks_seq[start-1:end]
                        print('>' + pks_id + '|' + pks_gca + '|' + pks_clade + '\n' + dom_seq)

prot_to_mech = {}
prot_to_clade = {}
with open('../../../Product_Template_Analysis/Reference_Protein_Table_from_Liu_et_al_2015.txt') as of:
    for line in of:
        line = line.strip()
        ls = line.split()
        prot = ls[3]
        clade = ls[1]
        mech = ls[2]
        prot_to_mech[prot] = mech
        prot_to_clade[prot] = clade

with open('../../../Product_Template_Analysis/synthaser_calls_on_pks.json') as osf:
    d = json.load(osf)
    for pks in d:
        pks_id = pks['header']
        pks_seq = pks['sequence']
        mech = prot_to_mech[pks_id]
        clade = prot_to_clade[pks_id]
        for dom in pks['domains']:
            dt = dom['type']
            if dt == 'PT':
                start = dom['start']
                end = dom['end']
                dom_seq = pks_seq[start-1:end]
                print('>Ref|' + pks_id + '|' + mech + '|' + clade + '\n' + dom_seq) 
