###############################################
## Make all of the qsub files for PAIRADISE
## Author: Levon Demirdjian
## Last Updated: 10/23/2019
## Description: This script creates all of the 
## qsub files needed for running the PAIRADISE
## pipeline
###############################################


## Inputs: 
## 1) working_dir: Directory where all of the data are located
## 2) pop_nam: Population name
## 3) pop_info: The 3 column .txt file containing the sample ID's and RNA-seq read lengths.

## Read in the user inputs
working_dir=$1
pop_name=$2
pop_info=$3

## Make a directory to contain all of the qsub files
mkdir ${working_dir}/qsub_files

## Make other required directories
mkdir ${pop_name}
mkdir ASEvents
mkdir results
mkdir plots

## Count the number of samples
x=`wc -l < ${working_dir}/${pop_info}`

############
## STEP 1 ##
############

## Create .qsub files for step 1 (personalize, map, and assign)
for i in $(seq 1 $x); do

    sample_name=`echo $(sed "${i}q;d" ${working_dir}/${pop_info}) | awk '{print $1}'`
    grep -m 1 -A 9 "############# ${sample_name}" ${working_dir}/${pop_name}_pairadise_pipeline.qsub > ${working_dir}/qsub_files/pairadise_step1_${sample_name}.qsub

done


############
## STEP 2 ##
############

## Create a .qsub file for step 2 (joint annotation)
grep -m 1 -A 1 "pairadise_annotate" ${working_dir}/${pop_name}_pairadise_pipeline.qsub > ${working_dir}/qsub_files/pairadise_step2_annotation.qsub


############
## STEP 3 ##
############

## Create .qsub files for step 3 (individual counting)
for i in $(seq 1 $x); do

    sample_name=`echo $(sed "${i}q;d" ${working_dir}/${pop_info}) | awk '{print $1}'`
    echo Done with $sample_name
    grep -m 1 -A 1 "pairadise_count -o $pop_name/$sample_name" ${working_dir}/${pop_name}_pairadise_pipeline.qsub > ${working_dir}/qsub_files/pairadise_step3_${sample_name}.qsub

done


############
## STEP 4 ##
############

## Create a .qsub file for step 4 (merge counts)
grep -m 1 -A 1 "pairadise_count --merge" ${working_dir}/${pop_name}_pairadise_pipeline.qsub > ${working_dir}/qsub_files/pairadise_step4_merge.qsub


############
## STEP 5 ##
############

## Create a .qsub file for step 5 (statistical modeling)
grep -m 1 -A 10 '###### Running statistical model' ${working_dir}/${pop_name}_pairadise_pipeline.qsub > ${working_dir}/qsub_files/pairadise_step5_stat.qsub
