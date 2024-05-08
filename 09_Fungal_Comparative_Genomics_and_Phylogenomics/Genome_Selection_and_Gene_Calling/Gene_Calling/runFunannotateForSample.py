import os
import sys

genome_file = os.path.abspath(sys.argv[1])
sample_name = sys.argv[2]
output_dir = os.path.abspath(sys.argv[3]) + '/'
cpus = int(sys.argv[4])

os.system('mkdir ' + output_dir)

# step 1: process
cleaned_gf = output_dir + sample_name + '.cleaned.fa'
clean_cmd = 'funannotate clean -i ' + genome_file + ' -o ' + cleaned_gf
os.system(clean_cmd)

sort_gf = output_dir + sample_name + '.cleaned.sorted.fa'
sort_cmd = 'funannotate sort -i ' + cleaned_gf + ' -o ' + sort_gf
os.system(sort_cmd)

mask_gf = output_dir + sample_name + '.cleaned.sorted.masked.fa'
mask_cmd = 'funannotate mask --cpus ' + str(cpus) + ' -i ' + sort_gf + ' -o ' + mask_gf
os.system(mask_cmd)

# step 2: predict genes
results_dir = output_dir + 'gene_prediction_results/'
os.system('mkdir ' + results_dir)
predict_cmd = 'funannotate predict --cpus ' + str(cpus) + ' -i ' + mask_gf + ' -o ' + results_dir + ' -s ' + sample_name
os.system(predict_cmd)
