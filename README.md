# Ancestral_BGC_Expansions_Study

Code and datasets for analyses featured in Salamzade, Kalan, and Currie, 2024.

*Note, names of files/paths in the repo or in scripts/files were occassionally changed to be more interpretable for viewers of this git repository.*

The file **Supplementary_Figures.pdf** contains some additional figures related to the study. They might not be properly displayed on GitHub so downloading might be needed. 

## Directory stucture:

* **01_BacterialPhylum_Genome_Selection**: Code for selection of species-level representative genomes from five bacterial phyla (Actinomycetota, Bacillota, Cyanobacteriota, Myxococcota, and Pseudomonadota) for investigation of BGC content and associations with other genomic/ecological factors (Table S1).
* **02_BacterialAndArchaealOrder_Genome_Selection**: Code for selection of genus-level representative genomes from orders for the bacterial and archaeal timetree analysis (Table S2).
* **03_BacterialPhylum_Phylogeny_Construction**: Code and example IQ-TREE commands for construction of phylum level phylogenies based on 74 conserved proteins across bacteria from GToTree.
* **04_BacterialAndArchaealOrder_TimeTree_Construction**: Code and scripts for constructing the timetree of bacteria and archae order representatives using MCMCTree.
* **05_BacterialAndArchaealOrder_Ornstein-Uhlenbeck_BGC_Count_Shift**: Code to run Ornstein-Uhlenbeck analysis to detect expansions in median BGC-ome size across bacterial order representatives.
* **06_Systematic_Bacterial_SpeciesLevel_Annotations**: Code and example commands/scripts for various trait annotations of species-level representative genomes (e.g. antiSMASH, CAZy/dbCAN, etc.).
---
* **07_Actinomycetota_Comparative_Genomics**: Code pertaining to primary Actinomycetota analysis on core set of 133 high-quality genomes for which OrthoFinder was run.
* **08_Actinomycetota_Comprehensive_MIBiG_ResFam_and_RifamycinRelated_Annotations**: Code pertaining to more large-scale Actnomycetota analyses on all ~33,000 genomes in the phylum in GTDB R214 - primarily around identifying MIBiG BGCs, Resfam resistance markers, rifamycin resistance, and fai-based detection of rifamycin BGC homologs. 
---
* **09_Fungal_Comparative_Genomics_and_Phylogenomics**: Code pertaining to fungal genomics analyses:
  * **Genome_Selection_and_Gene_Calling**: Code for genome selection of representative fungal genomes. Code for selecting gene calling information from NCBI databases or performing gene calling using funannotate/MetaEuk.
  * **Phylogenomics**: Code for comprehensive phylogeny construction of all 252 fungal genome representatives and 4 outgroup genomes.
  * **TimeTree_Construction**: Code for construction of timetree using MCMCTree and subset of representative genomes from the fungal kingdom.
  * **Annotations**: Code for performing annotation of BGCs using antiSMASH, CAZy/dBCAN homologs, and conidium regulators.
  * **Key_BGC_Gene_Typing_and_Ortholog_Grouping**: Code for typing and assessing conservation of key BGC genes (NRP, PK, and terpene biosynthesis) across fungal taxonomic clades based on synthaser and OrthoFinder analyses.
  * **BGCome_Size_Comparisons_Between_Taxonomic_Partitions**: Code for comparing counts of BGCs by antiSMASH between fungal taxonomic groups.
  * **Taxonomic_Clades**: A directory containing files for different taxonomic or morphological genome listings.
  * **Correlation_Analyses**: Code for correlation analyses of BGC-ome size with other factors, including genome size and number of distinct CAZy homologs.
  * **OrthoFinder_Comparative_Genomics_and_GWAS_Analyses**: Code pertaining to OrthoFinder analysis, comparative genomics between basal and BGC-enriched Pezizomycotina, and Dikarya GWAS for HOGs associated with increased NRPS+PKS-ome size.
