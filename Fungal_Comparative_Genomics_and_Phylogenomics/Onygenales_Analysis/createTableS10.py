import os
import sys

as_file = 'Ony_AntiSMASH_Stats.txt'
ov_file = 'Ony_Overview_File.txt'
meta_file = 'Ony_MetaInfo.txt'

meta = {}
with open(meta_file) as omf:
    for line in omf:
        line = line.strip()
        ls = line.split('\t')
        meta[ls[1].split('.')[0]] = ls[2:]

annot_method = {}
with open(ov_file) as ovf:
    for line in ovf:
        line = line.strip()
        ls = line.split('\t')
        annot_method[ls[0]] = [ls[1], ls[2]]

with open(as_file) as oaf:
    for i, line in enumerate(oaf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: 
            print('\t'.join([ls[0]] + ['Taxonomy', 'Annotation method', 'Morphology', 'Morphology source'] + [x.replace('_', ' ') for x in ls[1:]]))
        else:
            print('\t'.join([ls[0]] + annot_method[ls[0]] + meta[ls[0]] + ls[1:]))
