import os
import sys
from Bio import SeqIO
#implements the GARP:FIMNKY criterion for comp hetero per taxon suggested by Munoz-Gomez et al. 2018
from Bio import SeqIO, AlignIO
import sys, operator
import numpy as np
import pandas as pd
from scipy import stats
import statistics

def assessZ(seq_fasta):
    flat_ratios = []
    seqaln = SeqIO.index(seq_fasta, "fasta")
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
        flat_ratios.append(ratio)
    return(statistics.median(flat_ratios))

indir = 'Protein_Alignments_FurtherFiltered_GappySequencesRemoved/' # Directory of protein alignments post filtering of >10% gappy/ambiguous sites.
medians = []
for f in os.listdir(indir):
    if f.endswith('.faa'):
        medians.append(assessZ(indir + f))

print(statistics.median(medians))
