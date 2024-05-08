import os
import sys
from collections import defaultdict

sandpiper_file = 'Individual_Results/Median_RelAbundance_vs_HostProp.txt'
protcount_file = 'Individual_Results/Protein_Counts.txt'
antismash_stats_file = 'Individual_Results/Merged_Relaxed_and_Strict_Stats.txt'
species_file = 'Individual_Results/Species_Level_Selections.txt'
basic_stats_file = 'Individual_Results/Basic_Metrics.txt'
phispy_file = 'Individual_Results/PhiSpy_Summary.txt'
cazy_file = 'Individual_Results/CAZy_Counts_per_GCA.txt'
adherence_file = 'Individual_Results/Top_Hits_Blastp_Adherence.txt'
ko_file = 'Individual_Results/KOFam_Summary.txt'
myxo_ko_file = 'Individual_Results/Myxo_Social_and_Photosynthesis_Information.txt'
myxo_spo_file = 'Individual_Results/Myxo_Sporulation_Presence.txt'
actino_ssgb_file = 'Individual_Results/GCAs_with_SsgB_Based_on_HMM.txt'
actino_class_file = 'Individual_Results/Actinomycetia_GCAs.txt'
cyano_mc_file = 'Individual_Results/Cyano_Multicellularity_and_Sporulation_Presence.txt'
aa_freqs_file = 'Individual_Results/AA_Frequencies.txt'
transposon_file = 'Individual_Results/Transposon_Counts.txt'

ssgb = set([])
with open(actino_ssgb_file) as oasf:
    for line in oasf:
        line = line.strip()
        ssgb.add(line)

actinomycetia = set([])
with open(actino_class_file) as ocf:
    for line in ocf:
        line = line.strip()
        actinomycetia.add(line)

cyano_mc_genes = defaultdict(set)
cyano_spo_genes = defaultdict(set)
with open(cyano_mc_file) as ocmf:
    for line in ocmf:
        line = line.strip()
        ls = line.split('\t')
        if ls[-1] == 'Multicellularity':
            cyano_mc_genes[ls[0]].add(ls[1])
        elif ls[-1] == 'Sporulation':
            cyano_spo_genes[ls[0]].add(ls[1])
        elif ls[-1] == 'Sporulation/Multicellularity':
            cyano_spo_genes[ls[0]].add(ls[1])
            cyano_mc_genes[ls[0]].add(ls[1])

cazy_counts_total = defaultdict(int)
cazy_counts_distinct = defaultdict(int)
with open(cazy_file) as ocf:
    for line in ocf:
        line = line.strip()
        ls  = line.split('\t')
        gca = ls[0].split('.')[0]
        cazy_counts_total[gca] = int(ls[1])
        cazy_counts_distinct[gca] = int(ls[2])

adherence_counts = defaultdict(set)
with open(adherence_file) as oaf:
    for line in oaf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('.')[0]
        for hit in ls[2:]:
            hit_label = hit.split('(')[0]
            adherence_counts[gca].add(hit_label)

prot_counts = {}
with open(protcount_file) as opf:
    for line in opf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('.')[0]
        prot_counts[gca] = ls[1]

phispy_stats = {}
with open(phispy_file) as opf:
    for line in opf:
        line = line.strip()
        ls = line.split('\t')
        phispy_stats[ls[0]] = [float(ls[1]), float(ls[2])]

gca_to_sp = {}
gca_to_gen = {}
gca_to_phy = {}
with open(species_file) as osf:
    for line in osf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[-1].split('.')[0]
        gca_to_sp[gca] = ls[2]
        gca_to_gen[gca] = ls[1]
        gca_to_phy[gca] = ls[0]

gca_to_gc = {}
gca_to_thermal = {}
with open(basic_stats_file) as ogsf:
    for i, line in enumerate(ogsf):
        line = line.strip()
        ls = line.split('\t')
        if i == 0: continue
        gca = ls[0].split('.')[0]
        gca_to_gc[gca] = ls[3]
        gca_to_thermal[gca] = ls[-1]

oxy_counts = defaultdict(int)
with open(ko_file) as okf:
    for i, line in enumerate(okf):
        line = line.strip()
        if i == 0: continue
        ls = line.split('\t')
        if ls[-1] != 'Oxidative Phosphorylation': continue
        gca = ls[0]
        oxy_counts[gca] += 1

myxo_social_counts = defaultdict(int)
with open(myxo_ko_file) as okf:
    for i, line in enumerate(okf):
        line = line.strip()
        if i == 0: continue
        ls = line.split('\t')
        gca = ls[0]
        if ls[-1] != 'Photosynthesis':
            myxo_social_counts[gca] += 1

myxo_spo_genes = defaultdict(set)
with open(myxo_spo_file) as omsf:
    for line in omsf:
        line = line.strip()
        ls = line.split('\t')
        myxo_spo_genes[ls[0]].add(ls[1])

species_to_sandpiper_info = defaultdict(lambda: ['not-detected', '0.0', '0.0', '0', '0', '0'])
with open(sandpiper_file) as ospf:
    for i, line in enumerate(ospf):
        # phylum  species host_associated_category     median_relative_abundance       proportion_host_associated      total_metagenome_count
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        spinfo = ls[2:]
        sp = ls[1][3:]
        species_to_sandpiper_info[sp] = spinfo

aa_headers = []
aa_freqs = {}
with open(aa_freqs_file) as oaff:
    for i, line in enumerate(oaff):
        line = line.strip()
        ls = line.split('\t')
        if i == 0:
            aa_headers = []
            for j, val in enumerate(ls[1:]):
                if j == 12 or j == 18: continue
                aa_headers.append(val)
        else:
            gca = ls[0].split('.')[0]
            vals = []
            for j, val in enumerate(ls[1:]):
                if j == 12 or j == 18: continue
                vals.append(val)
            aa_freqs[gca] = vals 

transposon_count = defaultdict(int)
with open(transposon_file) as otf:
    for line in otf:
        line = line.strip()
        ls = line.split('\t')
        transposon_count[ls[0]] = int(ls[1])

print('\t'.join(['gca', 'relaxed_bgcome_genome_prop', 'relaxed_bgcome_size', 'strict_bgcome_size', 'strict_nrps_or_pks_ome_size', 'phylum', 'genus', 'species', 'host_associated_category', 'median_relative_abundance', 'proportion_host_associated', 'total_metagenome_count', 'host_metagenome_count', 'environment_metagenome_count', 'genome_size', 'protein_count', 'gc', 'ideal_temp', 'phage_sum', 'transposon_count', 'cazy_total_count', 'cazy_distinct_count', 'oxidative_phospho_count', 'adherence_count', 'Cyanobacteriota_multicellularity_genes', 'Cyanobacteriota_sporulation_genes', 'Actinomycetota_ssgB_status', 'Myxococcota_social_genes', 'Myxococcota_sporulation_genes'] + aa_headers))
with open(antismash_stats_file) as oasf:
    for i, line in enumerate(oasf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('.')[0]
        bgc_sum_relaxed, bgc_sum_strict, _, nrps_or_pks_sum_strict, _, _, genome_size, bgc_genome_prop_relaxed = ls[1:9]
        
        sp = gca_to_sp[gca]
        genus = gca_to_gen[gca]
        phylum = gca_to_phy[gca]

        if phylum in set(['Chloroflexota', 'Bacteroidota']): continue  # note should actually be 'Pseudomonadota'
        
        if phylum == 'Pseudomondota': phylum = 'Pseudomonadota'
        prot_count = float(prot_counts[gca])
        gca_phage_sum = str(phispy_stats[gca][1])
        gca_trans_count = str(transposon_count[gca])
        gca_oxy_count = str(oxy_counts[gca])#/prot_count)
        gca_caz_total_count = str(cazy_counts_total[gca])#/prot_count)
        gca_caz_distinct_count = str(cazy_counts_distinct[gca])#/prot_count)
        gca_adh_count = str(len(adherence_counts[gca]))#/prot_count)
        
        if genus in set(['GCA-2687015', 'JACZVU01', 'UBA4427']): continue

        gca_aa_freqs = aa_freqs[gca]
        prot_count = str(prot_count)
        top_mg_category, median_rel_abd, prop_host_associated, total_mg_count, host_mg_count, env_mg_count = species_to_sandpiper_info[sp]
        gc = gca_to_gc[gca]
        thermal_pref = gca_to_thermal[gca]
        actino_ssgb_status = 'NA'
        if gca in actinomycetia:
            if gca in ssgb:
                actino_ssgb_status = 'Actinomycetia with SsgB'
            else:
                actino_ssgb_status = 'Actinomycetia without SsgB'
        elif phylum == 'Actinomycetota':
            if gca in ssgb:
                actino_ssgb_status = 'Other Actinomycetota with SsgB'
            else:
                actino_ssgb_status = 'Other Actinomycetota without SsgB'
                
        myxo_soc = 'NA'
        myxo_spo = 'NA'
        if phylum == 'Myxococcota':
            myxo_soc = '0'
            myxo_spo = '0'
            if gca in myxo_social_counts:
                myxo_soc = str(myxo_social_counts[gca])
            if gca in myxo_spo_genes:
                myxo_spo = str(len(myxo_spo_genes[gca]))

        cyano_mc = 'NA'
        cyano_spo = 'NA'
        if phylum == 'Cyanobacteriota':
            cyano_spo = '0'
            cyano_mc = '0'
            if gca in cyano_mc_genes:
                cyano_mc = str(len(cyano_mc_genes[gca]))
            if gca in cyano_spo_genes:
                cyano_spo = str(len(cyano_spo_genes[gca]))

        print('\t'.join([gca, bgc_genome_prop_relaxed, bgc_sum_relaxed, bgc_sum_strict, nrps_or_pks_sum_strict, phylum, genus, sp, top_mg_category, median_rel_abd, prop_host_associated, total_mg_count, host_mg_count, env_mg_count, genome_size, prot_count, gc, thermal_pref, gca_phage_sum, gca_trans_count, gca_caz_total_count, gca_caz_distinct_count, gca_oxy_count, gca_adh_count, cyano_mc, cyano_spo, actino_ssgb_status,  myxo_soc, myxo_spo] + aa_freqs[gca])) 
