# Ancestral_BGC_Expansions_Study

Code and select datasets for analyses featured in Salamzade, Kalan and Currie, 2025.

*Note, names of files/paths in the repo or in scripts/files were occassionally changed to be more interpretable for viewers of this git repository.*

The file **Supplementary_Figures.pdf** contains some additional figures related to the study. They might not be properly displayed on GitHub so downloading might be needed. 

Additional datasets related to this project can be found on Zenodo at: https://zenodo.org/uploads/15168076

## Directory stucture:

* **01_BacterialAndArchaealOrder_Genome_Selection**: Code for selection of genus-level representative genomes from orders for the bacterial and archaeal timetree analysis (Supplementary Table S1).
* **02_BacterialPhylum_Genome_Selection**: Code for selection of species-level representative genomes from five bacterial phyla (Actinomycetota, Bacillota, Cyanobacteriota, Myxococcota, and Pseudomonadota) for investigation of BGC content and associations with other genomic/ecological factors (Supplementary Table S3).
* **03_BacterialPhylum_Phylogeny_Construction**: Code and example IQ-TREE commands for construction of phylum level phylogenies based on 74 conserved proteins across bacteria from GToTree.
* **04_BacterialAndArchaealOrder_TimeTree_Construction**: Code and scripts for constructing the timetree of bacteria and archae order representatives using MCMCTree.
* **05_BacterialAndArchaealOrder_Ornstein-Uhlenbeck_BGC_Count_Shift**: Code to run Ornstein-Uhlenbeck analysis to detect expansions in median BGC-ome size across bacterial order representatives.
* **06_Systematic_Bacterial_SpeciesLevel_Annotations_and_Associations**: Code and example commands/scripts for various trait annotations of species-level representative genomes (e.g. antiSMASH, CAZy/dbCAN, etc.).
---
* **07_Actinomycetota_Comparative_Genomics**: Code pertaining to primary Actinomycetota analysis on core set of 133 high-quality genomes for which OrthoFinder was run.
* **08_Actinomycetota_Comprehensive_MIBiG_ResFam_and_RifamycinRelated_Annotations**: Code pertaining to more large-scale Actnomycetota analyses on all ~33k genomes in the phylum in GTDB R214 - primarily around identifying MIBiG BGCs, Resfam resistance markers, rifamycin resistance, and fai-based detection of rifamycin BGC homologs. 
---
* **09_Fungal_Comparative_Genomics_and_Phylogenomics**: Code pertaining to fungal genomics analyses:
  * **Genome_Selection_and_Gene_Calling**: Code for genome selection of representative fungal genomes. Code for selecting gene calling information from NCBI databases or performing gene calling using funannotate (Supplementary Table S7).
  * **Phylogenomics**: Code for comprehensive phylogeny construction of all 312 fungal genome representatives + 4 outgroup genomes.
  * **TimeTree_Construction**: Code for construction of timetree using MCMCTree and subset of representative genomes from the fungal kingdom.
  * **Annotations**: Code for performing annotation of BGCs using antiSMASH & synthaser, CAZy/dBCAN profile HMMs, and HI domains.
  * **PKS_and_NRPS_from_Basal_in_Context_of_Characterized_MIBIG**: Code/data used to perform NaPDoS2-based analysis of whether basal/early-diverging fungal genomes harbor novel PKS/NRPS sequences to what has been characterized and is on MIBiG.
  * **BGCome_Size_Comparisons_Between_Taxonomic_Partitions**: Code for comparing counts of BGCs by antiSMASH between select fungal taxonomic groups.
  * **Taxonomic_Clades**: A directory containing files for different taxonomic or morphological genome listings.
  * **Correlation_Analyses**: Code for correlation analyses of BGC-ome size with other factors, including genome size and number of distinct CAZy homologs.
  * **OrthoFinder_Comparative_Genomics_and_GWAS_Analyses**: Code/select-data pertaining to OrthoFinder analysis, comparative genomics between BGC-enriched and BGC-poor dikaryotic clades, and fungal kingdom-wide GWAS for HOGs associated with increased NRPS+PKS-ome size.
  * **HI_Domain_Architectures_and_Phylogenies**: Code and select data pertaining to deeper dives into HI domain containing proteins.
  * **HI_Enrichment_Analysis**: Code for performing assessment of whether BGCs have elevated rates of HI domain containing proteins in BGC-rich Pezizomycotina genomes.
  * **Pangenome_Expansion_Analysis**: Code/data pertaining to use of psaps for assessing genomic fluidity of fungal clades. The psaps code can be found at: https://github.com/raufs/psaps 
  * **Aspergillus_Investigations**: Select data for phylogenomics, annotations, and psaps analysis of *Aspergillus* genomes with regards to HI protein associations with BGCs and genome fluidity.
  * **Ideograms**: Code for creating ideograms showcasing the location of BGC regions and HI domain containing proteins along select scaffolds from two dikaryotic genomes.
* **10_Recent_Horizontal_Gene_Transfer_Assessment**: Code and select data for assessing signatures of recent horizontal transfer in bacteria, archaeae, and fungi.

