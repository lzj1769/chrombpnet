#!/bin/bash
##TODO flank_size and ref_fasta files in scripts 

cell_line=HEPG2
data_type="ATAC"
neg_shift=4

date=$(date +'%m.%d.%Y')
#setting=naked_bias_4_$neg_shift"_shifted_"$data_type"_"$date
cur_file_name="hepg2_atac_script_new.sh"
setting=naked_bias_4_4_shifted_ATAC_07.29.2021
### SIGNAL INPUT

naked_bam=/oak/stanford/groups/akundaje/projects/enzymatic_bias_correction/pipeline_out/atac/SRR072187/call-filter/shard-0/execution/SRR072187.4_1.merged.nodup.no_chrM_MT.bam
in_bam=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/merged_data/HEPG2.atac.filt.merged.bam
overlap_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/25b3429e-5864-4e8d-a475-a92df8938887/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz
#.gz file?
idr_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/25b3429e-5864-4e8d-a475-a92df8938887/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
is_filtered=True
samtools_flag=None

blacklist_region=$PWD/data/all_three_blacklists.bed
chrom_sizes=$PWD/data/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/$cell_line
data_dir=$PWD/$cell_line/data
naked_data_dir=$PWD/naked_dna_scripts/data
output_dir=$PWD/$cell_line/$setting


### MODEL PARAMS

gpu=0
filters=500 
n_dil_layers=8
seed=1234 
model_name=model 
neg_dir=$main_dir/neg_data
flank_size=1057

## CREATE DIRS

if [[ -d $main_dir ]] ; then
    echo "main director already exists"
else
    mkdir $main_dir
fi

if [[ -d $output_dir ]] ; then
    echo "output director already exists"
else
    mkdir $output_dir
fi

###  BIAS  INPUTS

bias_json=""
bias_weights=""
neg_bed_train=$neg_dir"/bpnet.inputs.train.flank300.0.negatives.bed"
neg_bed_test=$neg_dir"/bpnet.inputs.test.0.negatives.bed"

step2_bias_json=$output_dir/bias_fit_on_signal_step2/model.0.arch
step2_bias_weights=$output_dir/bias_fit_on_signal_step2/model.0.weights


### CREATE BIGWIGS ON INVIVO DATA

if [[ -d $data_dir ]] ; then
    echo "skipping bigwig creation"
else
    mkdir $data_dir
    ./main_scripts/preprocess.sh $in_bam $data_dir $samtools_flag $is_filtered $data_type $neg_shift
    cp $PWD/naked_dna_scripts/$cur_file_name $data_dir
fi


### CREATE FOLDER INVIVO TILEDB AND RUN TILEDB

if [[ -d $data_dir/tiledb ]] ; then
    echo "skipping tiledb"
else
    mkdir $data_dir/tiledb
    echo -e "dataset\tnegatives_peak\tidr_peak\toverlap_peak\tambig_peak\tcount_bigwig_unstranded_5p\n"$cell_line"\t"$neg_dir/bpnet.inputs.all.negatives.bed"\t"$idr_peak"\t"$overlap_peak"\t"$blacklist_region"\t"$data_dir/shifted.sorted.bam.bpnet.unstranded.bw > $data_dir/tiledb/inputs.tsv
    echo -e "overlap_peak\tbed_summit_from_last_col\nnegatives_peak\tbed_summit_from_last_col\nidr_peak\tbed_summit_from_last_col\nambig_peak\tbed_no_summit\ncount_bigwig_unstranded_5p\tbigwig" > $data_dir/tiledb/attribs.txt
    ./main_scripts/db_ingest.sh  $data_dir/tiledb/inputs.tsv $data_dir/tiledb/db $chrom_sizes $data_dir/tiledb/attribs.txt
    cp $PWD/naked_dna_scripts/$cur_file_name $data_dir/tiledb
fi

### CREATE BIGWIGS ON NAKED DATA

if [[ -d $naked_data_dir ]] ; then
    echo "skipping naked dna bigwig creation"
else
    mkdir $naked_data_dir
    ./main_scripts/preprocess.sh $naked_bam $naked_data_dir $samtools_flag $is_filtered $data_type $neg_shift
    cp $PWD/naked_dna_scripts/$cur_file_name $naked_data_dir
fi


### CREATE FOLDER NAKED TILEDB AND RUN TILEDB

if [[ -d $naked_data_dir/tiledb ]] ; then
    echo "skipping naked tiledb"
else
    mkdir $naked_data_dir/tiledb
    echo -e "dataset\tidr_peak\toverlap_peak\tambig_peak\tcount_bigwig_unstranded_5p\n"$cell_line"\t"$idr_peak"\t"$overlap_peak"\t"$blacklist_region"\t"$naked_data_dir/shifted.sorted.bam.bpnet.unstranded.bw > $naked_data_dir/tiledb/inputs.tsv
    echo -e "overlap_peak\tbed_summit_from_last_col\nidr_peak\tbed_summit_from_last_col\nambig_peak\tbed_no_summit\ncount_bigwig_unstranded_5p\tbigwig" > $naked_data_dir/tiledb/attribs.txt
    ./main_scripts/db_ingest.sh  $naked_data_dir/tiledb/inputs.tsv $naked_data_dir/tiledb/db $chrom_sizes $naked_data_dir/tiledb/attribs.txt
    cp $PWD/naked_dna_scripts/$cur_file_name $naked_data_dir/tiledb
fi


### STEP 1 -  TRAIN BIAS MODEL


if [[ -d $naked_data_dir/naked_bias_model_step1 ]] ; then
    echo "skipping step 1 - directory present "
else

    if test -z "$bias_json" 
    then
    	mkdir $naked_data_dir/naked_bias_model_step1
        grep -w -v "chr1" $overlap_peak > $naked_data_dir/naked_bias_model_step1/train.bed
        grep -w "chr1" $overlap_peak > $naked_data_dir/naked_bias_model_step1/test.bed

        bash main_scripts/get_loss_weights.sh $naked_data_dir/tiledb/db "chr10" "overlap_peak" "count_bigwig_unstranded_5p" $cell_line $flank_size $naked_data_dir/naked_bias_model_step1/counts_loss_weight.txt
        counts_loss_weight_step1=`cat $naked_data_dir/naked_bias_model_step1/counts_loss_weight.txt`
 
   	echo -e "counts_loss_weight\t"$counts_loss_weight_step1"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $naked_data_dir/naked_bias_model_step1/params.txt
    	params=$naked_data_dir/naked_bias_model_step1/params.txt
    	for fold in 0
    	do
            ./main_scripts/naked_bias_model_step1/train.sh $fold $gpu $model_name $seed $naked_data_dir/naked_bias_model_step1 $params  $naked_data_dir/tiledb/db $cell_line profile_bpnet_dnase $naked_data_dir/naked_bias_model_step1/train.bed
            ./main_scripts/naked_bias_model_step1/predict.sh $fold $gpu $model_name $seed  $naked_data_dir/naked_bias_model_step1  $naked_data_dir/tiledb/db $cell_line $chrom_sizes $naked_data_dir/naked_bias_model_step1/test.bed
            ./main_scripts/naked_bias_model_step1/score.sh $naked_data_dir/naked_bias_model_step1 $model_name $fold $cell_line $seed
    	done
    	cp $PWD/naked_dna_scripts/$cur_file_name $naked_data_dir/naked_bias_model_step1
        bias_json=$naked_data_dir/naked_bias_model_step1/model.0.arch
        bias_weights=$naked_data_dir/naked_bias_model_step1/model.0.weights
    else
        echo "skipping step1 - input bias model given"
    fi

fi


bias_json=$naked_data_dir/naked_bias_model_step1/model.0.arch
bias_weights=$naked_data_dir/naked_bias_model_step1/model.0.weights

### STEP 2 - FIT BIAS MODEL ON SIGNAL


if [[ -d $output_dir/bias_fit_on_signal_step2 ]] ; then
    echo "skipping step 2"
else
    mkdir $output_dir/bias_fit_on_signal_step2

    bash main_scripts/get_loss_weights.sh $data_dir/tiledb/db "chr10" "overlap_peak" "count_bigwig_unstranded_5p" $cell_line $flank_size $output_dir/bias_fit_on_signal_step2/counts_loss_weight.txt
    counts_loss_weight_step2=`cat $output_dir/bias_fit_on_signal_step2/counts_loss_weight.txt`
    counts_loss_weight_step3=$counts_loss_weight_step2

    echo -e "json_string\t"$bias_json"\nweights\t"$bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step2"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/bias_fit_on_signal_step2/params.txt
    params=$output_dir/bias_fit_on_signal_step2/params.txt
    for fold in 0
    do
        ./main_scripts/bias_fit_on_signal_step2/train.sh $fold $gpu $model_name $seed $output_dir/bias_fit_on_signal_step2 $params  $data_dir/tiledb/db $cell_line $PWD/naked_dna_scripts/main_scripts/bias_fit_on_signal_step2/signal_from_bias.py
        ./main_scripts/bias_fit_on_signal_step2/predict.sh $fold $gpu $model_name $seed  $output_dir/bias_fit_on_signal_step2  $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/bias_fit_on_signal_step2/score.sh $output_dir/bias_fit_on_signal_step2 $model_name $fold $cell_line $seed
    done
    cp $PWD/naked_dna_scripts/$cur_file_name $output_dir/bias_fit_on_signal_step2
fi

counts_loss_weight_step2=`cat $output_dir/bias_fit_on_signal_step2/counts_loss_weight.txt`
counts_loss_weight_step3=$counts_loss_weight_step2



### STEP 3 - FIT BIAS AND SIGNAL MODEL

if [[ -d $output_dir/final_model_step3_new ]] ; then
    echo "skipping step 3"
else
    mkdir $output_dir/final_model_step3_new
    echo -e "json_string\t"$step2_bias_json"\nweights\t"$step2_bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step3"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/final_model_step3_new/params.txt
    params=$output_dir/final_model_step3_new/params.txt
    for fold in 0
    do
        ./main_scripts/final_model_step3_new/train.sh $fold $gpu $model_name $seed $output_dir/final_model_step3_new $params  $data_dir/tiledb/db $cell_line $PWD/main_scripts/final_model_step3_new/profile_bpnet_dnase_with_bias.py
        ./main_scripts/final_model_step3_new/predict.sh $fold $gpu $model_name $seed  $output_dir/final_model_step3_new  $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/final_model_step3_new/score.sh $output_dir/final_model_step3_new $model_name $fold $cell_line $seed
    done
    cp $PWD/naked_dna_scripts/$cur_file_name $output_dir/final_model_step3_new
fi


## UNPLUG MODEL

if [[ -d $output_dir/final_model_step3_new/unplug ]] ; then
    echo "skipping unplugging"
else
    mkdir $output_dir/final_model_step3_new/unplug
    unplug_bias_json=$output_dir/final_model_step3_new/model.0.arch
    unplug_bias_weights=$output_dir/final_model_step3_new/model.0.weights
    echo -e "json_string\t"$unplug_bias_json"\nweights\t"$unplug_bias_weights"\ncounts_loss_weight\t"$counts_loss_weight_step3"\nprofile_loss_weight\t1\nfilters\t"$filters"\nn_dil_layers\t"$n_dil_layers > $output_dir/final_model_step3_new/unplug/params.txt

    params=$output_dir/final_model_step3_new/unplug/params.txt

    for fold in 0
    do
        CUDA_VISIBLE_DEVICES=$gpu python ./main_scripts/unplug_new/get_model_with_bias_unplugged.py --model_params $params --outf $output_dir/final_model_step3_new/unplug/$model_name.$fold.hdf5 
        ./main_scripts/unplug/predict.sh $fold $gpu $model_name $seed $output_dir/final_model_step3_new/unplug $data_dir/tiledb/db $cell_line $chrom_sizes
        ./main_scripts/unplug/score.sh $output_dir/final_model_step3_new/unplug $model_name $fold $cell_line $seed
    done
    cp $PWD/naked_dna_scripts/$cur_file_name $output_dir/final_model_step3_new/unplug
fi


### GET INTERPRETATIONS

if [[ -d $data_dir/$cell_line"_idr_split" ]] ; then
    echo "skipping creating idr splits for interpretation"
else
    mkdir  $data_dir/$cell_line"_idr_split" 
    zcat $idr_peak | shuf  > $data_dir/$cell_line"_idr_split/temp.txt"
    split -l 10000 $data_dir/$cell_line"_idr_split/temp.txt" $data_dir/$cell_line"_idr_split/x"
    rm  $data_dir/$cell_line"_idr_split/temp.txt"
fi



if [[ -d $output_dir/final_model_step3_new/unplug/deepshap ]] ; then
    echo "skipping interpretations"
else
    mkdir $output_dir/final_model_step3_new/unplug/deepshap
    bed_file=$data_dir/$cell_line"_idr_split"

    for fold in 0
    do
        ./main_scripts/interpret/interpret.sh $output_dir/final_model_step3_new/unplug/$model_name.$fold.hdf5 $bed_file xaa $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3_new/unplug/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret.sh $output_dir/final_model_step3_new/unplug/$model_name.$fold.hdf5 $bed_file xab $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3_new/unplug/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret.sh $output_dir/final_model_step3_new/unplug/$model_name.$fold.hdf5 $bed_file xac $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3_new/unplug/deepshap $cell_line $gpu $fold
    done

    python $PWD/main_scripts/interpret/combine_shap_pickle.py --source $output_dir/final_model_step3_new/unplug/deepshap --target $output_dir/final_model_step3_new/unplug/deepshap --type 20k
    cp $PWD/naked_dna_scripts/$cur_file_name $output_dir/final_model_step3_new/unplug/deepshap

fi

fold=0
bed_file=$data_dir/$cell_line"_idr_split"
./main_scripts/interpret/interpret.sh $output_dir/final_model_step3_new/unplug/$model_name.$fold.hdf5 $bed_file xac $data_dir/tiledb/db $chrom_sizes $output_dir/final_model_step3_new/unplug/deepshap $cell_line $gpu $fold
python $PWD/main_scripts/interpret/combine_shap_pickle.py --source $output_dir/final_model_step3_new/unplug/deepshap --target $output_dir/final_model_step3_new/unplug/deepshap --type 20k

modisco_sig_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/SIGNAL/
if [[ -d $modisco_sig_dir/$cell_line ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_sig_dir/$cell_line
fi

modisco_sig_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/SIGNAL/
if [[ -d $modisco_sig_dir/$cell_line/$setting"_new"/ ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_sig_dir/$cell_line/$setting"_new"/
    modisco_dir_final=$modisco_sig_dir/$cell_line/$setting"_new"/
    cp  $cell_line/$setting/final_model_step3_new/unplug/deepshap/20K.fold0.deepSHAP $modisco_dir_final
fi

#naked model interpretations

if [[ -d $naked_data_dir/naked_bias_model_step1/deepshap ]] ; then
    echo "skipping interpretations"
else
    mkdir $naked_data_dir/naked_bias_model_step1/deepshap
    bed_file=$data_dir/$cell_line"_idr_split"

    for fold in 0
    do
        ./main_scripts/interpret/interpret_weight.sh $naked_data_dir/naked_bias_model_step1/$model_name.$fold $bed_file xaa $naked_data_dir/tiledb/db $chrom_sizes $naked_data_dir/naked_bias_model_step1/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret_weight.sh $naked_data_dir/naked_bias_model_step1/$model_name.$fold $bed_file xab $naked_data_dir/tiledb/db $chrom_sizes $naked_data_dir/naked_bias_model_step1/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret_weight.sh $naked_data_dir/naked_bias_model_step1/$model_name.$fold $bed_file xac $naked_data_dir/tiledb/db $chrom_sizes $naked_data_dir/naked_bias_model_step1/deepshap $cell_line $gpu $fold
    done

    python $PWD/main_scripts/interpret/combine_shap_pickle.py --source $naked_data_dir/naked_bias_model_step1/deepshap --target $naked_data_dir/naked_bias_model_step1/deepshap --type 20k
    cp $PWD/naked_dna_scripts/$cur_file_name $naked_data_dir/naked_bias_model_step1/deepshap

fi



modisco_sig_dir=/oak/stanford/groups/akundaje/projects/chrombpnet_paper/importance_scores/BIAS/
if [[ -d $modisco_sig_dir/naked_dna ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_sig_dir/naked_dna
fi


if [[ -d $modisco_sig_dir/naked_dna/$setting/ ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_sig_dir/naked_dna/$setting/
    modisco_dir_final=$modisco_sig_dir/naked_dna/$setting/
    cp  $naked_data_dir/naked_bias_model_step1/deepshap/20K.fold0.deepSHAP $modisco_dir_final
fi



# naked model after tuning interpretations

if [[ -d $output_dir/bias_fit_on_signal_step2/deepshap ]] ; then
    echo "skipping interpretations"
else
    mkdir $output_dir/bias_fit_on_signal_step2/deepshap
    bed_file=$data_dir/$cell_line"_idr_split"

    for fold in 0
    do
        ./main_scripts/interpret/interpret_weight.sh $output_dir/bias_fit_on_signal_step2/$model_name.$fold $bed_file xaa $data_dir/tiledb/db $chrom_sizes $output_dir/bias_fit_on_signal_step2/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret_weight.sh $output_dir/bias_fit_on_signal_step2/$model_name.$fold $bed_file xab $data_dir/tiledb/db $chrom_sizes $output_dir/bias_fit_on_signal_step2/deepshap $cell_line $gpu $fold
        ./main_scripts/interpret/interpret_weight.sh $output_dir/bias_fit_on_signal_step2/$model_name.$fold $bed_file xac $data_dir/tiledb/db $chrom_sizes $output_dir/bias_fit_on_signal_step2/deepshap $cell_line $gpu $fold
    done

    python $PWD/main_scripts/interpret/combine_shap_pickle.py --source $output_dir/bias_fit_on_signal_step2/deepshap --target $output_dir/bias_fit_on_signal_step2/deepshap --type 20k
    cp $PWD/naked_dna_scripts/$cur_file_name $output_dir/bias_fit_on_signal_step2/deepshap

fi


if [[ -d $modisco_sig_dir/$cell_line/$setting"_biasfit"/ ]] ; then
    echo "modisco dir already exists"
else
    mkdir $modisco_sig_dir/$cell_line/$setting"_biasfit"/
    modisco_dir_final=$modisco_sig_dir/$cell_line/$setting"_biasfit"/
    cp  $output_dir/bias_fit_on_signal_step2/deepshap/20K.fold0.deepSHAP $modisco_dir_final
fi

### RUN MODISCO


