import os 
import sys
from Bio import SeqIO
from collections import defaultdict

ancestral_conserved_file = 'Outgroup_Conserved_Positions_and_Alleles.txt'
simple_actino_families_dir = 'Other_Actinomycetia/'
part2_file = 'Mycobacteriaceae.txt' # Rotate: 'Micromonosporaceae.txt' # 'Mycobacteriaceae.txt' #'Micromonosporaceae.txt' #'Pseudonocardiaceae.txt'
part1_file = 'Streptomycetaceae.txt' 
core_faa_file = 'Core_Alignment_Concatenated.faa'

ancestral = {}
with open(ancestral_conserved_file) as oacf:
    for line in oacf:
        line = line.strip()
        ls = line.split('\t')
        ancestral[int(ls[0])] = ls[1]

part1_genomes = set([])
with open(part1_file) as op1f:
    for line in op1f:
        line = line.strip()
        part1_genomes.add(line)

part2_genomes = set([])
with open(part2_file) as op2f:
    for line in op2f:
        line = line.strip()
        part2_genomes.add(line)

print('\t'.join(['Family', 'Focal', 'Part1', 'Part2', 'Stat', 'Total_Count', 'AABA_and_BABA_Count', 'ABAA_and_BBAA_Count']))
for f in os.listdir(simple_actino_families_dir):
    focal_fam_file = simple_actino_families_dir + f
    focal_genomes = set([])
    with open(focal_fam_file) as offf:
        for line in offf:
            line = line.strip()
            focal_genomes.add(line)

    seqs = []
    index_to_genome = {}
    with open(core_faa_file) as ocff:
        for i, rec in enumerate(SeqIO.parse(ocff, 'fasta')):
            seqs.append(list(str(rec.seq)))
            index_to_genome[i] = rec.id

    # hypothesis for phylogeny: aaba > abaa & bbaa > baba
    focal_bbaa_or_aaba = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    focal_baba_or_abaa = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for pos, als in enumerate(zip(*seqs)):
        if pos in ancestral:
            anc_al = ancestral[pos]
            als = list(als)
            focal_genomes_alleles = {}
            part1_genomes_alleles = {}
            part2_genomes_alleles = {}
            for i, al in enumerate(als):
                genome = index_to_genome[i]
                if genome in focal_genomes:
                    focal_genomes_alleles[genome] = al.upper()
                elif genome in part1_genomes:
                    part1_genomes_alleles[genome] = al.upper()
                elif genome in part2_genomes:
                    part2_genomes_alleles[genome] = al.upper()

            for fg in focal_genomes_alleles:
                for p1 in part1_genomes_alleles:
                    for p2 in part2_genomes_alleles:
                        if focal_genomes_alleles[fg] != anc_al and part1_genomes_alleles[p1] != anc_al and focal_genomes_alleles[fg] == part1_genomes_alleles[p1] and part2_genomes_alleles[p2] == anc_al and focal_genomes_alleles[fg] != '-':
                            focal_bbaa_or_aaba[fg][p1][p2] += 1 
                        elif focal_genomes_alleles[fg] == anc_al and part1_genomes_alleles[p1] == anc_al and part2_genomes_alleles[p2] != anc_al and part2_genomes_alleles[p2] != '-':
                            focal_bbaa_or_aaba[fg][p1][p2] += 1
                        elif focal_genomes_alleles[fg] != anc_al and part1_genomes_alleles[p1] == anc_al and focal_genomes_alleles[fg] == part2_genomes_alleles[p2] and part2_genomes_alleles[p2] != anc_al and focal_genomes_alleles[fg] != '-': 
                            focal_baba_or_abaa[fg][p1][p2] += 1
                        elif focal_genomes_alleles[fg] == anc_al and part1_genomes_alleles[p1] != anc_al and part2_genomes_alleles[p2] == anc_al and part1_genomes_alleles[p1] != '-':
                            focal_baba_or_abaa[fg][p1][p2] += 1

    fam = f.split('.txt')[0]
    for fg in focal_genomes:
        for p1 in part1_genomes:
            for p2 in part2_genomes:
                bbaa_or_aaba = focal_bbaa_or_aaba[fg][p1][p2]
                baba_or_abaa = focal_baba_or_abaa[fg][p1][p2]
                tot = bbaa_or_aaba + baba_or_abaa
                stat = 'NA'
                if tot > 0:
                    stat = (bbaa_or_aaba - baba_or_abaa)/tot
                print('\t'.join([str(x) for x in [fam, fg, p1, p2, stat, tot, bbaa_or_aaba, baba_or_abaa]]))
