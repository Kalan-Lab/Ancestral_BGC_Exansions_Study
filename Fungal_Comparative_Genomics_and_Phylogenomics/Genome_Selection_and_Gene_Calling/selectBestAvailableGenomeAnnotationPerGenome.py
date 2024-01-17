import os
import sys

all_genomes_dir = os.path.abspath('All_Genomes/') + '/'
ncbi_genbank_dir = os.path.abspath('Genome_GFFs_in_NCBI/') + '/'
ncbi_refseq_dir = os.path.abspath('Genome_GFFs_in_NCBI_RefSeq/') + '/'
ncbi_refseq_fasta_dir = os.path.abspath('Genome_Fastas_in_NCBI_RefSeq/') + '/'
funannotate_dir = os.path.abspath('/researchdrive_files/Personal/Rauf/Multicellularity_and_BGCs/Fungal_Funannotate_Results_Dir/') + '/'
prob_gcas_file = 'Problematic_NCBI_GCAs_or_GCFs.txt' # these are GFF files which antiSMASH v7 had difficulties processing

prob_gcas = set([])
with open(prob_gcas_file) as opgf:
    for line in opgf:
        line = line.strip()
        prob_gcas.add(line)

ncbi_gbk_available = {}
funannotate_annot_available = set([])
ncbi_refseq_fasta = {}
ncbi_gffs = {} 
for f in os.listdir(ncbi_refseq_dir):
    gca_with_version = '_'.join(f.split('_')[:2])
    gca = gca_with_version.split('.')[0].replace('GCF_', 'GCA_')
    ncbi_gbk_available[gca] = gca_with_version.replace('GCA_', 'GCF_')
    ncbi_gffs[gca] = ncbi_refseq_dir + f

for f in os.listdir(ncbi_refseq_fasta_dir):
    gca_with_version = '_'.join(f.split('_')[:2])
    gca = gca_with_version.split('.')[0].replace('GCF_', 'GCA_')
    assert(ncbi_gbk_available[gca] == gca_with_version)
    ncbi_refseq_fasta[gca] = ncbi_refseq_fasta_dir + f

for f in os.listdir(ncbi_genbank_dir):
    gca_with_version = '_'.join(f.split('_')[:2])
    gca = gca_with_version.split('.')[0]
    if gca in ncbi_gbk_available: continue
    ncbi_gbk_available[gca] = gca_with_version
    ncbi_gffs[gca] = ncbi_genbank_dir + f

fun_gffs = {}
fun_fastas = {}
for gca in os.listdir(funannotate_dir):
    gff_file = funannotate_dir + gca + '/gene_prediction_results/predict_results/' + gca + '.gff3'
    fasta_file = funannotate_dir + gca + '/' + gca + '.cleaned.sorted.masked.fa'
    if os.path.isfile(gff_file) and os.path.getsize(gff_file) > 100:
        funannotate_annot_available.add(gca)
        fun_gffs[gca] = gff_file 
        assert(os.path.isfile(fasta_file))
        fun_fastas[gca] = fasta_file

for f in os.listdir(all_genomes_dir):
    gca_with_version = '_'.join(f.split('_')[:2])
    gca = gca_with_version.split('.')[0]
    if gca == 'GCA_001013415': continue
    ncbi_gbk_avail = 'False'
    ncbi_gbk_full = None
    if gca in ncbi_gbk_available: 
        ncbi_gbk_avail = 'True'
        ncbi_gbk_full = ncbi_gbk_available[gca]
    
    if ncbi_gbk_full != None:
        gca_with_version = ncbi_gbk_full
        #assert(ncbi_gbk_full == gca_with_version)

    fun_avail = 'False'
    if gca in funannotate_annot_available: 
        fun_avail = 'True'
    
    fasta_file = all_genomes_dir + f
    if gca in ncbi_refseq_fasta:
        fasta_file = ncbi_refseq_fasta[gca]

    gff_file = None
    annotation_method = None
    if gca in ncbi_gffs and not gca in prob_gcas:
        gff_file = ncbi_gffs[gca]
        if 'GCF_' in gca_with_version: annotation_method ='NCBI RefSeq'
        else: annotation_method= 'NCBI GenBank'
    elif gca in fun_gffs:
        annotation_method = 'funannotate'
        gff_file = fun_gffs[gca]
        fasta_file = fun_fastas[gca]
    elif gca == 'GCA_000412615':
        annotation_method = 'MetaEuk'
        gff_file = "/researchdrive_files/Personal/Rauf/Multicellularity_and_BGCs/James_Genomes_MetaEuk_GFFs_Reformatted/GCA_000412615.1_RESULT_.gff"
    elif gca == 'GCA_017911155':
        annotation_method = 'MetaEuk'
        gff_file = '/researchdrive_files/Personal/Rauf/Multicellularity_and_BGCs/James_Genomes_MetaEuk_GFFs_Reformatted/GCA_017911155.1_RESULT_.gff'
    else:
        print('ISSUE: %s' % gca)
        sys.exit(1)
    print('\t'.join([gca, annotation_method, gca_with_version, ncbi_gbk_avail, fun_avail, fasta_file, gff_file]))
