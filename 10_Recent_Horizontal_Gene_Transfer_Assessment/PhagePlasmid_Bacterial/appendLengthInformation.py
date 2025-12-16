import os
import sys

bgc_length = {}
with open('BGC_Information.txt') as obif:
    for line in obif:
        line = line.strip()
        ls = line.split('\t')
        bgc_length[ls[1]] = ls[5]

ignore = set(['Chloroflexota', 'Bacteroidota'])
with open('MGE_Overlap.txt') as omof:
    for line in omof:
        line = line.strip()
        ls = line.split("\t")
        if ls[2] in ignore: continue 
        print('\t'.join(ls + [bgc_length[ls[1]]]))
