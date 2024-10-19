import os
import sys
import logging
import traceback
import subprocess
import numpy
from Bio import SeqIO
import gzip
from collections import defaultdict

def n50_calc(genome_file):
	#Solution adapted from dinovski:
	#https://gist.github.com/dinovski/2bcdcc770d5388c6fcc8a656e5dbe53c
	lengths = []
	seq = ""	
	if genome_file.endswith('.gz'):
		with gzip.open(genome_file, 'rt') as fasta:
			for line in fasta:
				if line.startswith('>'):
					if seq != "":
						lengths.append(len(seq))
					seq = ""
				else:
					seq += line.strip()
		if seq != "":
			lengths.append(len(seq))
	else:	
		with open(genome_file) as fasta:
			for line in fasta:
				if line.startswith('>'):
					if seq != "":
						lengths.append(len(seq))
					seq = ""
				else:
					seq += line.strip()
		if seq != "":
			lengths.append(len(seq))

	## sort contigs longest>shortest
	all_len=sorted(lengths, reverse=True)
	csum=numpy.cumsum(all_len)

	n2=int(sum(lengths)/2)

	# get index for cumsum >= N/2
	csumn2=min(csum[csum >= n2])
	ind=numpy.where(csum == csumn2)
	n50 = all_len[int(ind[0])]
	return(n50)


taxid_to_name = {}
taxid_to_genus = {}
with open('TaxonKit_Names.txt') as otn:
    for line in otn:
        line = line.strip()
        ls = line.split('\t')
        taxid_to_name[ls[0]] = ls[1]
        species_index = None
        for i, val in enumerate(ls[2].split(';')):
            if val == 'species':
                species_index = i
        species = ''
        for i, val in enumerate(ls[1].split(';')):
            if i == species_index:
                species = val
        genus = species.replace('[', '').replace(']', '').replace("'", '').split()[0]
        assert(not '_' in genus)
        taxid_to_genus[ls[0]] = genus
        
gca_to_genus = {}
gca_to_taxonomy = {}
with open('All_UFCG_Representative_Selections.with_NCBI_TaxIDs.txt') as of:
    for line in of:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[1].split('.')[0]
        gca_to_taxonomy[gca] = taxid_to_name[ls[2]]
        gca_to_genus[gca] = taxid_to_genus[ls[2]]

with open('Previous_Overview_File.NCBI_TaxIDs.txt') as of:
    for line in of:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('.')[0].replace('GCF_', 'GCA_')
        gca_to_taxonomy[gca] = taxid_to_name[ls[1]]
        gca_to_genus[gca] = taxid_to_genus[ls[1]]

accounted_for = set([])
print('\t'.join(['gca', 'name', 'taxonomy', 'annotation_method', 'version_info', 'N50', 'fasta_file', 'gff_file']))
with open('Modified_Previous_Overview_File.txt') as ompof:
    for i, line in enumerate(ompof):
        if i == 0: continue
        line = line.strip()
        gca, taxonomy, annot_method, version_info, ncbi_annot_avail, funannotate_annot_avail, fasta_file, gff_file = line.split('\t')
        accounted_for.add(gca)
        
        gff_file = gff_file.replace('/home/salamzade/Multicellularity_and_BGC_Abundances/Polished_Fungi/Genome_GFFs_in_NCBI_RefSeq/', '/workspace/lab/kalanlab/salamzar/multicellularity_and_bgcs/calcaneus/Multicellularity_and_BGC_Abundances/Polished_Fungi/Genome_GFFs_in_NCBI_RefSeq/')
        fasta_file = fasta_file.replace('/home/salamzade/Multicellularity_and_BGC_Abundances/Polished_Fungi/Genome_Fastas_in_NCBI_RefSeq/', '/workspace/lab/kalanlab/salamzar/multicellularity_and_bgcs/calcaneus/Multicellularity_and_BGC_Abundances/Polished_Fungi/Genome_Fastas_in_NCBI_RefSeq/')
        
        gff_file = gff_file.replace('/home/salamzade/Multicellularity_and_BGC_Abundances/Polished_Fungi/Genome_GFFs_in_NCBI/', '/workspace/lab/kalanlab/salamzar/multicellularity_and_bgcs/calcaneus/Multicellularity_and_BGC_Abundances/Polished_Fungi/Genome_GFFs_in_NCBI/')
        fasta_file = fasta_file.replace('/home/salamzade/Multicellularity_and_BGC_Abundances/Polished_Fungi/All_Genomes/', '/workspace/lab/kalanlab/salamzar/multicellularity_and_bgcs/calcaneus/Multicellularity_and_BGC_Abundances/Polished_Fungi/All_Genomes_LinksCorrected/')

        gff_file = gff_file.replace('/researchdrive_files/Personal/Rauf/Multicellularity_and_BGCs/Fungal_Funannotate_Results_Dir/', '/workspace/lab/kalanlab/salamzar/multicellularity_and_bgcs/researchdrive/Multicellularity_and_BGCs/Fungal_Funannotate_Results_Dir/')
        fasta_file = fasta_file.replace('/researchdrive_files/Personal/Rauf/Multicellularity_and_BGCs/Fungal_Funannotate_Results_Dir/', '/workspace/lab/kalanlab/salamzar/multicellularity_and_bgcs/researchdrive/Multicellularity_and_BGCs/Fungal_Funannotate_Results_Dir/')
        fasta_file = fasta_file.replace('/All_Genomes/', '/All_Genomes_LinksCorrected/')

        taxonomy = gca_to_taxonomy[gca]
        name = gca_to_genus[gca] + '_' + version_info
        n50 = n50_calc(fasta_file)
        new_ls = [gca, name, taxonomy, annot_method, version_info, str(n50), fasta_file, gff_file]
        print('\t'.join(new_ls))

gca_fna_dir = os.path.abspath('../UFCG_Selections/Genomes/') + '/'
gcf_fna_dir = os.path.abspath('../Gene_Calling_for_New_Genomes/GCF_FASTAs/') + '/'
gca_gff_dir = os.path.abspath('../Gene_Calling_for_New_Genomes/GCA_GFFs/') + '/'
gcf_gff_dir = os.path.abspath('../Gene_Calling_for_New_Genomes/GCF_GFFs/') + '/'
fun_dir = os.path.abspath('/workspace/local/rauf/multicellularity_and_bgcs/funannotate_new_genomes/final_results/') + '/'

fun_fna = {}
fun_gff = {}
for subdir, dirs, files in os.walk(fun_dir):
    for file in files:
        full_file = os.path.join(subdir, file)
        if full_file.endswith('.cleaned.sorted.masked.fa'):
            gca = full_file.split('/')[-1].split('.cleaned.sorted.masked.fa')[0]
            fun_fna[gca] = full_file 
        elif 'predict_results/' in full_file and full_file.endswith('.gff3'):
            gca = full_file.split('/')[-1].split('.')[0]
            fun_gff[gca] = full_file
        
gca_fna = {}
gca_gff = {}
for f in os.listdir(gca_fna_dir):
    gca_fna[f.split('.')[0]] = gca_fna_dir + f

for f in os.listdir(gca_gff_dir):
    gca_gff[f.split('.')[0]] = gca_gff_dir + f

gcf_fna = {}
gcf_gff = {}
for f in os.listdir(gcf_fna_dir):
    gcf_fna[f.split('.')[0].replace('GCF_', 'GCA_')] = gcf_fna_dir + f

for f in os.listdir(gcf_gff_dir):
    gcf_gff[f.split('.')[0].replace('GCF_', 'GCA_')] = gcf_gff_dir + f

as_gcfs = set([])
for f in os.listdir('AntiSMASH_Results_GCF/'):
    gca = f
    as_gcfs.add(gca)

prob_gca_gffs = set(['GCA_000524645', 'GCA_000827475', 'GCA_000218685', 'GCA_000442785'])

with open('GCAs.txt') as ogf:
    for line in ogf:
        version = line.strip()
        gca = line.strip().split('.')[0]
        gff_file = "Not found"
        fasta_file = "Not found"
        annot_method = "NA"
        if gca in as_gcfs:
            annot_method = "NCBI RefSeq"
            version = version.replace('GCA_', 'GCF_')
            fasta_file = gcf_fna[gca]
            gff_file = gcf_gff[gca]
        elif gca in gca_gff and not gca in prob_gca_gffs:
            annot_method = 'NCBI GenBank'
            fasta_file = gca_fna[gca]
            gff_file = gca_gff[gca]
        else:
            annot_method = 'funannotate'
            try:
                fasta_file = fun_fna[gca]
                gff_file = fun_gff[gca]
            except:
                sys.stderr.write('ISSUES!\n%s\n--------\n'% gca)
                continue

        taxonomy = gca_to_taxonomy[gca]
        name = gca_to_genus[gca] + '_' + version

        n50 = "NA"
        if fasta_file != 'Not found':
            n50 = n50_calc(fasta_file)
        
        new_ls = [gca, name, taxonomy, annot_method, version, str(n50), fasta_file, gff_file]
        print('\t'.join(new_ls))

