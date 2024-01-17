# implements the GARP:FIMNKY criterion for comp hetero per taxon suggested by Munoz-Gomez et al. 2018
# code adapted from: https://github.com/Tancata/phylo/blob/master/rank_taxa_by_aa_compo.py

from Bio import SeqIO, AlignIO
import sys, operator
import numpy as np
import pandas as pd
from scipy import stats

#for this dataset, ratio of about 1 seems to be a reasonable cutoff between low and high GC.
highGC = []
lowGC = []
ratios = {}
flat_ratios = []
flat_zeds = []
col_zed = {}

seqaln = SeqIO.index(sys.argv[1], "fasta") # individual protein alignment
ratio_midpoint = float(sys.argv[2]) # ratio midpoint determined per phylum for distinguishing GC rich and poor organisms
for taxon in seqaln:
    garp = 0
    fimnky = 0
    sequence = str(seqaln[taxon].seq)
    for char in sequence:
        if char in 'GARP':
            garp += 1
        elif char in 'FIMNKY':
            fimnky += 1
    ratio = float(garp)/float(fimnky)
    if ratio > ratio_midpoint:
        highGC.append(taxon)
    else:
        lowGC.append(taxon)
    ratios[taxon] = ratio
    flat_ratios.append(ratio)

#now compute Z for each column in the alignment

alignment = AlignIO.read(sys.argv[1], "fasta")
for col in range(len(str(alignment[0].seq))):
    GARP_high = 0.0
    GARP_low = 0.0
    FIMNKY_high = 0.0
    FIMNKY_low = 0.0
    highGC_taxa = float(len(highGC))
    lowGC_taxa = float(len(lowGC))
    for rec in alignment:
        char = str(rec.seq)[col]
        if rec.id in highGC:
            if char in 'GARP':
                GARP_high += 1.0/highGC_taxa
            elif char in 'FIMNKY':
                FIMNKY_high += 1.0/highGC_taxa
        elif rec.id in lowGC:
            if char in 'GARP':
                GARP_low += 1.0/lowGC_taxa
            elif char in 'FIMNKY':
                FIMNKY_low += 1.0/lowGC_taxa
        else:
            print("Taxon " + str(rec.id) + " couldn't be assigned to high or low GC...")
            quit()

    zed = (FIMNKY_low - FIMNKY_high) + (GARP_high - GARP_low)
    flat_zeds.append(zed) 
    col_zed[col] = zed

#now write alignments in which top X% of sites by zed have been removed - say 50% for now
sorted_zeds = sorted(col_zed.items(), key=operator.itemgetter(1))
num_cols = len(sorted_zeds)
cols_to_take = int(float(num_cols)/2.0)
selected_cols = []
for i in range(cols_to_take):
    selected_cols.append(sorted_zeds[i][0])

for rec in seqaln:
    seq_to_print = ''
    for col in sorted(selected_cols):
        seq_to_print += str(seqaln[rec].seq)[col]
    print(">" + str(rec) + "\n" + str(seq_to_print) + "\n")
