import os
import sys
from collections import defaultdict

sandpiper_file = '../../General_CAZy_to_BGC_Analysis/Median_RelAbundance_vs_HostProp.txt'
protcount_file = '../../AntiSMASH_For_All_Species_Reps/Protein_Counts.txt'
antismash_stats_file = '../../AntiSMASH_For_All_Species_Reps/Merged_Relaxed_and_Strict_Stats.txt'
species_file = '../../AntiSMASH_For_All_Species_Reps/Species_Level_Selections.txt'
basic_stats_file = '../../AntiSMASH_For_All_Species_Reps/Basic_Metrics.txt'
genomad_file = '../../AntiSMASH_For_All_Species_Reps/geNomad_Summary.txt'
cazy_file = '../../AntiSMASH_For_All_Species_Reps/dbCAN_Analysis/HMM_Based/CAZy_Counts_per_GCA.txt'
adherence_file = '../../AntiSMASH_For_All_Species_Reps/Adherence_Search/Top_Hits_Blastp.txt'
ko_file = '../../AntiSMASH_For_All_Species_Reps/KOFam_Summary.txt'
myxo_ko_file = '../../Myxo/KO_Annotations/Myxo_Social_and_Photosynthesis_Information.txt'
myxo_spo_file = '../../Myxo/Sporulation_Annotation/Myxo_Sporulation_Presence.txt'
actino_ssgb_file = '../../Actino/SsgB_Homolog_Finding/GCAs_with_SsgB_Based_on_HMM.txt'
actino_class_file = '../../Actino/SsgB_Homolog_Finding/Actinomycetia_GCAs.txt'
firmi_spor_file = '../../Firmi/Sporulation/GCA_COG_Counts.txt'
firmi_multicel_file = '../../Firmi/Sporulation/GCA_MC_COG_Counts.txt'
cyano_mc_file = '../../Cyano/Multicellularity_and_Sporulation_Annotations/Cyano_Multicellularity_and_Sporulation_Presence.txt'
aa_freqs_file = '../../AntiSMASH_For_All_Species_Reps/AA_Frequencies.txt'
transposon_file = '../../AntiSMASH_For_All_Species_Reps/ISFinder_Analysis/Transposon_Counts.txt'

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

bac_spo_counts = {}
with open(firmi_spor_file) as ofsf:
    for line in ofsf:
        line = line.strip()
        ls = line.split('\t')
        bac_spo_counts[ls[0].split('.')[0]] = ls[1]

bac_mc_counts = {}
with open(firmi_multicel_file) as ofmf:
    for line in ofmf:
        line = line.strip()
        ls = line.split('\t')
        bac_mc_counts[ls[0].split('.')[0]] = ls[1]

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

genomad_stats = {}
with open(genomad_file) as opf:
    for i, line in enumerate(opf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        genomad_stats[ls[0]] = [float(ls[1]), float(ls[2]), float(ls[3])]

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

myxo_photos_counts = defaultdict(int)
myxo_social_counts = defaultdict(int)
with open(myxo_ko_file) as okf:
    for i, line in enumerate(okf):
        line = line.strip()
        if i == 0: continue
        ls = line.split('\t')
        gca = ls[0]
        if ls[-1] == 'Photosynthesis':
            myxo_photos_counts[gca] += 1
        else:
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

print('\t'.join(['gca', 'relaxed_bgcome_genome_prop', 'relaxed_bgcome_size', 'strict_bgcome_size', 'strict_nrps_or_pks_ome_size', 'phylum', 'genus', 'species', 'host_associated_category', 'median_relative_abundance', 'proportion_host_associated', 'total_metagenome_count', 'host_metagenome_count', 'environment_metagenome_count', 'genome_size', 'protein_count', 'gc', 'ideal_temp', 'phage_plasmid_sum', 'transposon_count', 'cazy_total_count', 'cazy_distinct_count', 'oxidative_phospho_count', 'adherence_count', 'Cyanobacteriota_multicellularity_genes', 'Cyanobacteriota_sporulation_genes', 'Actinomycetota_ssgB_status', 'Myxococcota_social_genes', 'Myxococcota_sporulation_genes'] + aa_headers))
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

        if phylum in set(['Chloroflexota', 'Bacteroidota']): continue 

        # note should actually be 'Pseudomonadota'
        if phylum == 'Pseudomondota': phylum = 'Pseudomonadota'
        prot_count = float(prot_counts[gca])
        gca_phage_plasmid_sum = str(genomad_stats[gca][2])
        gca_trans_count = str(transposon_count[gca])
        gca_oxy_count = str(oxy_counts[gca])#/prot_count)
        gca_caz_total_count = str(cazy_counts_total[gca])#/prot_count)
        gca_caz_distinct_count = str(cazy_counts_distinct[gca])#/prot_count)
        gca_adh_count = str(len(adherence_counts[gca]))#/prot_count)
        
        photo_status = 'Not Photosynthetic'
        if phylum == 'Myxococcota':
            if gca in myxo_photos_counts:
                pc = myxo_photos_counts[gca]
                if pc >= 15:
                    photo_status = 'Photosynthetic'
        if phylum == 'Cyanobacteriota':
           photo_status = 'Photosynthetic'

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

        bac_spo = 'NA'
        bac_mc = 'NA'
        if phylum == 'Bacillota':
            bac_mc = '0'
            ac_spo = '0'
            if gca in bac_mc_counts:
                bac_mc = bac_mc_counts[gca]
            if gca in bac_spo_counts:
                bac_spo = bac_spo_counts[gca]
            
        cyano_mc = 'NA'
        cyano_spo = 'NA'
        if phylum == 'Cyanobacteriota':
            cyano_spo = '0'
            cyano_mc = '0'
            if gca in cyano_mc_genes:
                cyano_mc = str(len(cyano_mc_genes[gca]))
            if gca in cyano_spo_genes:
                cyano_spo = str(len(cyano_spo_genes[gca]))

        print('\t'.join([gca, bgc_genome_prop_relaxed, bgc_sum_relaxed, bgc_sum_strict, nrps_or_pks_sum_strict, phylum, genus, sp, top_mg_category, median_rel_abd, prop_host_associated, total_mg_count, host_mg_count, env_mg_count, genome_size, prot_count, gc, thermal_pref, gca_phage_plasmid_sum, gca_trans_count, gca_caz_total_count, gca_caz_distinct_count, gca_oxy_count, gca_adh_count, cyano_mc, cyano_spo, actino_ssgb_status,  myxo_soc, myxo_spo] + aa_freqs[gca])) 
