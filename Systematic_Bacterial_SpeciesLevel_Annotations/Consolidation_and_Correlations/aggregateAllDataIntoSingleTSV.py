import os
import sys
from collections import defaultdict

sandpiper_file = 'Median_RelAbundance_vs_HostProp.txt'
protcount_file = 'Protein_Counts.txt'
antismash_stats_file = 'AntiSMASH_For_All_Species_Reps/AntiSMASH_Stats.Updated.txt'
species_file = 'Species_Level_Selections.txt'
basic_stats_file = 'Basic_Metrics.txt'
phispy_file = 'PhiSpy_Summary.txt'
cazy_file = 'Top_Hits_Blastp.txt'
adherence_file = 'Top_Hits_Blastp.txt'
ko_file = 'KOFam_Summary.txt'
myxo_ko_file = 'Myxo_Social_and_Photosynthesis_Information.txt'
actino_ssgb_file = 'GCAs_with_SsgB_Based_on_HMM.txt'
actino_class_file = 'SsgB_Homolog_Finding/Actinomycetia_GCAs.txt'
firmi_spor_file = 'Bacillota_GCA_COG_Counts.txt'
firmi_multicel_file = 'Bacillota_GCA_MC_COG_Counts.txt'
cyano_mc_file = 'Cyano_Multicellularity_Data.txt'
aa_freqs_file = 'AA_Frequencies.txt'
transposon_file = 'Transposon_Counts.txt'

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

spor = {}
with open(firmi_spor_file) as ofsf:
    for line in ofsf:
        line = line.strip()
        ls = line.split('\t')
        spor[ls[0].split('.')[0]] = ls[1]

multicel = {}
with open(firmi_multicel_file) as ofmf:
    for line in ofmf:
        line = line.strip()
        ls = line.split('\t')
        multicel[ls[0].split('.')[0]] = ls[1]

cyano_mc_counts = {}
with open(cyano_mc_file) as ocmf:
    for line in ocmf:
        line = line.strip()
        ls = line.split('\t')
        cyano_mc_counts[ls[0]] = ls[1]

cazy_counts = defaultdict(set)
with open(cazy_file) as ocf:
    for line in ocf:
        line = line.strip()
        ls  = line.split('\t')
        gca = ls[0].split('.')[0]
        for hit in ls[2:]:
            hit_label = hit.split('|')[-1]
            cazy_counts[gca].add(hit_label)

adherence_counts = defaultdict(set)
with open(adherence_file) as oaf:
    for line in oaf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('.')[0]
        for hit in ls[2:]:
            hit_label = hit.split('(')[0]
            adherence_counts[gca].add(hit_label)

genome_sizes = {}
with open(antismash_stats_file) as oasf:
    for i, line in enumerate(oasf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        aid, phylum, genus, bgc, genome_size, bgc_genome_prop, _, nrps_genome_prop, pks_genome_prop, nrps_or_pks_genome_prop = ls
        gca = aid.split('.')[0]
        genome_sizes[gca] = genome_size

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
        if not ls[0] in genome_sizes: continue
        phispy_stats[ls[0]] = [float(ls[1]), float(ls[2])]

gca_to_sp = {}
with open(species_file) as osf:
    for line in osf:
        line = line.strip()
        ls = line.split('\t')
        gca = ls[-1].split('.')[0]
        gca_to_sp[gca] = ls[2]

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

photos_counts = defaultdict(int)
social_counts = defaultdict(int)
with open(myxo_ko_file) as okf:
    for i, line in enumerate(okf):
        line = line.strip()
        if i == 0: continue
        ls = line.split('\t')
        gca = ls[0]
        if ls[-1] == 'Photosynthesis':
            photos_counts[gca] += 1
        else:
            social_counts[gca] += 1

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

print('\t'.join(['gca', 'bgc_prop', 'bgc_sum', 'phylum', 'genus', 'species', 'host_associated_category', 'median_relative_abundance', 'proportion_host_associated', 'total_metagenome_count', 'host_metagenome_count', 'environment_metagenome_count', 'genome_size', 'protein_count', 'gc', 'ideal_temp', 'phage_sum', 'transposon_count', 'cazy_count', 'oxidative_phospho_count', 'photosynthesis_ability', 'adherence_count', 'Cyanobacteriota_multicellularity_genes', 'Actinomycetota_ssgb_status', 'Bacillota_sporulation_genes', 'Bacillot_multicellularity_genes', 'Myxococcota_social_genes'] + aa_headers))
with open(antismash_stats_file) as oasf:
    for i, line in enumerate(oasf):
        if i == 0: continue
        line = line.strip()
        ls = line.split('\t')
        gca = ls[0].split('.')[0]
        if not gca in genome_sizes: continue
        aid, phylum, genus, bgc, genome_size, bgc_genome_prop, _, nrps_genome_prop, pks_genome_prop, nrps_or_pks_genome_prop = ls 
        if phylum == 'Chloroflexota' or phylum == 'Bacteroidota': continue
        if phylum == 'Pseudomondota': phylum = 'Pseudomonadota'
        sp = gca_to_sp[gca]
        prot_count = float(prot_counts[gca])
        gca_phage_sum = str(phispy_stats[gca][1])
        gca_trans_count = str(transposon_count[gca])
        gca_oxy_count = str(oxy_counts[gca])#/prot_count)
        gca_caz_count = str(len(cazy_counts[gca]))#/prot_count)
        gca_adh_count = str(len(adherence_counts[gca]))#/prot_count)
        
        photo_status = 'Not Photosynthetic'
        if gca in photos_counts:
            pc = photos_counts[gca]
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
                
        bgc_sum = str(float(bgc_genome_prop)*float(genome_size))
        myxo_social_genes = 'NA'
        if phylum == 'Myxococcota':
            myxo_social_genes = '0'
            if gca in social_counts:
                myxo_social_genes = str(social_counts[gca])
        firmi_sporulation_genes = 'NA'
        if phylum == 'Bacillota':
            firmi_sporulation_genes = '0'
            if gca in spor:
                firmi_sporulation_genes = str(spor[gca])
        firmi_multicellularity_genes = 'NA'
        if phylum == 'Bacillota':
            if gca in multicel:
                firmi_multicellularity_genes = str(multicel[gca])
        
        cyano_mc = 'NA'
        if phylum == 'Cyanobacteriota':
            if gca in cyano_mc_counts:
                cyano_mc = cyano_mc_counts[gca]
            else:
                cyano_mc = '0'

        print('\t'.join([gca, bgc_genome_prop, bgc_sum, phylum, genus, sp, top_mg_category, median_rel_abd, prop_host_associated, total_mg_count, host_mg_count, env_mg_count, genome_size, prot_count, gc, thermal_pref, gca_phage_sum, gca_trans_count, gca_caz_count, gca_oxy_count, photo_status, gca_adh_count, cyano_mc, actino_ssgb_status, firmi_sporulation_genes, firmi_multicellularity_genes, myxo_social_genes] + aa_freqs[gca])) 
