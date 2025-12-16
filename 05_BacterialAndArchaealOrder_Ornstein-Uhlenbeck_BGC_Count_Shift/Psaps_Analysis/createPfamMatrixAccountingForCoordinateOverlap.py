import os
import sys
from collections import defaultdict
from ete3 import Tree
from Bio import SeqIO
from operator import itemgetter

faa_file = 'BGC_BiosynRuleBased_Proteins.faa'
pfam_annot_file = 'BGC_BiosynRuleBased_Protein_Pfam_Annotations_with_Coords.txt'
t = Tree("Renamed_Phylogeny.tre")

in_phylo = set([])
for n in t.traverse('postorder'):
    if n.is_leaf():
        in_phylo.add(n.name)

bit_arrays = {}
with open(faa_file) as off:
    for rec in SeqIO.parse(off, 'fasta'):
        bit_arrays[rec.id] = [0]*len(str(rec.seq))

all_pfams = set([])
all_samples = set([])
sample_pfams = defaultdict(set)


protein_dom_hits = defaultdict(list)
with open(pfam_annot_file) as odaf:
    for line in odaf:
        line = line.strip()
        prot, dom, start, end, score, evalue  = line.split('\t')
        protein_dom_hits[prot].append([dom, int(start), int(end), float(score)])

protein_dom_hits_retained = defaultdict(list)
for prot in protein_dom_hits:
    gca = prot.split('|')[0]
    coords_hit_already = set([])
    for di in sorted(protein_dom_hits[prot], key=itemgetter(3), reverse=True):
        overlapping_better_hits = 0
        dom_coords = set([])
        for pp in range(di[1], di[2]+1):
            dom_coords.add(pp)
            if pp in coords_hit_already:
                overlapping_better_hits += 1
        prop_overlapping = overlapping_better_hits/len(dom_coords)
        if prop_overlapping < 0.1:
            protein_dom_hits_retained[prot].append(di)
            coords_hit_already = coords_hit_already.union(dom_coords)
            
            all_pfams.add(di[0])
            all_samples.add(gca)
            sample_pfams[gca].add(di[0])

print('Pfam\t' + '\t'.join(sorted(all_pfams)))
for s in sorted(all_samples):
    if not s in in_phylo: continue
    row = [s]
    for pf in sorted(all_pfams):
        if pf in sample_pfams[s]:
            row.append('1')
        else:
            row.append('0')
    print('\t'.join(row))
