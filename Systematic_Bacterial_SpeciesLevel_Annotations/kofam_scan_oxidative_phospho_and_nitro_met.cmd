# in the following KOFamScan_DB_Filtered - only KOs belonging to the pathway Oxidative Phosphorylation or Nitrogen Metabolism were included for computational effiency.

exec_annotation -p KOFamScan_DB_Filtered/select_ko_profiles/ -k KOFamScan_DB_Filtered/ko_list --cpu 1 --tmp-dir prepTG_Results_21k/KOFamScan_Select_Tmp_Dir/GCA_905193035.1_ERR1600539_mag_bin.13.faa/ -o prepTG_Results_21k/KOFamScan_Select_Results/GCA_905193035.1_ERR1600539_mag_bin.13.faa.txt prepTG_Results_21k/Proteomes/GCA_905193035.1_ERR1600539_mag_bin.13.faa
