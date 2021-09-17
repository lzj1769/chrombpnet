flank=500
input_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/SIGNAL/
output_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_$flank/SIGNAL/
tomtom_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/tomtom/flank_len_$flank/SIGNAL/

model_name=4_4_shifted_ATAC_09.12.2021_bias_filters_500_new

#bash run_tomtom.sh SURAG $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh GM12878 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh HEPG2 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh IMR90 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh H1 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh K562 $model_name $input_dir $output_dir $tomtom_dir


input_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/BIAS/
output_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_$flank/BIAS/
tomtom_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/tomtom/flank_len_$flank/BIAS/

model_name=4_4_shifted_ATAC_09.12.2021_bias_filters_500

#bash run_tomtom.sh naked_dna $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh GM12878 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh HEPG2 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh IMR90 $model_name $input_dir $output_dir $tomtom_dir
#bash run_tomtom.sh H1 $model_name $input_dir $output_dir $tomtom_dir
bash run_tomtom.sh K562 $model_name $input_dir $output_dir $tomtom_dir
