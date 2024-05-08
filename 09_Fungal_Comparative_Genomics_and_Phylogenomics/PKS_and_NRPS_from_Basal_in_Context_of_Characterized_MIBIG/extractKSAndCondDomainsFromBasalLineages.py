import os
import sys
from collections import defaultdict
from Bio import SeqIO
import json 

all_pks_prots_file = 'synthaser_input_files/All_PKS_Full_Proteins.faa'
all_nrps_prots_file = 'synthaser_input_files/All_NRPS_Full_Proteins.faa'

cout = open('Condensation_Domain_Sequences.faa', 'w')
ksout = open('KetoacylSynthase_Domain_Sequences.faa', 'w')

genus_reps_file = 'Genus_Representative_GCAs.txt'

coi_set = set(['PKS|Type I|Non-reducing', 'PKS|Type I|Highly-reducing', 'PKS|Type I|Partially-reducing', 'Hybrid NRPS-PKS', 'Hybrid PKS-NRPS', 'NRPS'])
nrps_and_pks_prots_listing_file = 'NRPS_and_PKS_Proteins.txt'

pks_contam = set(['278', '279', '280', '281'])
nrps_contam = set(['144'])

redundant = set([])
redundant_names = set([])
with open(nrps_and_pks_prots_listing_file) as onaplf:
    for line in onaplf:
        line = line.strip()
        redundant.add(line)
        redundant_names.add(' '.join(line.split()[1:]))

genus_reps = set([])
with open(genus_reps_file) as ogrf:
    for line in ogrf:
        line = line.strip()
        genus_reps.add(line)

comp_clade_file = 'Clades/Dikarya_Zygomycota_Olpidium_and_Outgroups.txt'

complement = set([])
with open(comp_clade_file) as occf:
    for line in occf:
        line = line.strip()
        gca = '_'.join(line.split('_')[-2:])
        complement.add(gca)

id_to_gca = {}
redundant_pks_ids = set([])
with open(all_pks_prots_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        if rec.id in pks_contam: continue
        gca = rec.description.split('|')[1]
        id_to_gca[rec.id] = rec.description.split('|')[1]
        if ' '.join(rec.description.split()[1:]) in redundant_names:
            redundant_pks_ids.add(rec.id)
            
synth_files = ['synthaser_result_files/synthaser_calls_on_pks_p1.json', 'synthaser_result_files/synthaser_calls_on_pks_p2.json', 'synthaser_result_files/synthaser_calls_on_pks_p3.json', 'synthaser_result_files/synthaser_calls_on_pks_p4.json']

classifications = set([])
clade_subtype_gcas = defaultdict(lambda: defaultdict(set))
all_sts = set([])
for sf in synth_files:
    with open(sf) as osf:
        d = json.load(osf)
        for pks in d:
            pks_id = pks['header']
            if pks_id in pks_contam: continue
            pks_gca = id_to_gca[pks_id]
            if pks_gca in complement: continue
            if not pks_gca in genus_reps: continue
            if len(pks['classification']) > 0:
                pks_class = '|'.join(pks['classification'])
                if pks_class in coi_set:
                    for dom_id, dom in enumerate(pks['domains']):
                        if dom['type'] == 'C':
                            dom_start = dom['start']
                            dom_end = dom['end']
                            dom_seq = pks['sequence'][dom_start-0:dom_end]
                            cout.write('>' + pks_gca + '|' + pks_id + "|Domain_ID|" + str(dom_id) + "|Classification|" + pks_class + "|Domain|" + dom['domain'] + "\n" + dom_seq + "\n")
                        elif dom['type'] == 'KS':
                            dom_start = dom['start']
                            dom_end = dom['end']
                            dom_seq = pks['sequence'][dom_start-0:dom_end]
                            ksout.write('>' + pks_gca + '|' + pks_id + "|Domain_ID|" + str(dom_id) + "|Classification|" + pks_class + "|Domain|" + dom['domain'] + "\n" + dom_seq + "\n")

id_to_gca = {}
redundant_ids = set([])
with open(all_nrps_prots_file) as oapf:
    for rec in SeqIO.parse(oapf, 'fasta'):
        if rec.id in nrps_contam: continue
        gca = rec.description.split('|')[1]
        if rec.description in redundant:
            redundant_ids.add(rec.id)
        id_to_gca[rec.id] = rec.description.split('|')[1]

synth_files = ['synthaser_result_files/NRPS_synthaser_p1.json', 'synthaser_result_files/NRPS_synthaser_p2.json', 'synthaser_result_files/NRPS_synthaser_p3.json']

all_sts = set([])
for sf in synth_files:
    with open(sf) as osf:
        d = json.load(osf)
        for nrps in d:
            nrps_id = nrps['header']
            if nrps_id in nrps_contam: continue
            nrps_gca = id_to_gca[nrps_id]
            if nrps_gca in complement: continue
            if not nrps_gca in genus_reps: continue
            if len(nrps['classification']) > 0:
                nrps_class = '|'.join(nrps['classification'])
                if nrps_class in coi_set: 
                    for dom_id, dom in enumerate(nrps['domains']):
                        if dom['type'] == 'C':
                            dom_start = dom['start']
                            dom_end = dom['end']
                            dom_seq = nrps['sequence'][dom_start-0:dom_end]
                            cout.write('>' + nrps_gca + '|' + nrps_id + "|Domain_ID|" + str(dom_id) + "|Classification|" + nrps_class + "|Domain|" + dom['domain'] + "\n" + dom_seq + "\n")
                        elif dom['type'] == 'KS':
                            dom_start = dom['start']
                            dom_end = dom['end']
                            dom_seq = nrps['sequence'][dom_start-0:dom_end]
                            ksout.write('>' + nrps_gca + '|' + nrps_id + "|Domain_ID|" + str(dom_id) + "|Classification|" + nrps_class + "|Domain|" + dom['domain'] + "\n" + dom_seq + "\n")

cout.close()
ksout.close()
