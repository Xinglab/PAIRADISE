#!/bin/bash

#$ -cwd
#$ -o Output/array_out.$JOB_ID.stdout
#$ -e Output/array_out.$JOB_ID.stderr

. /etc/profile
module load python/2.7.14

## Pull the SGE_TASK_ID'th row of the parameter file as inputs
cur_line=$(sed "${SGE_TASK_ID}q;d" data/running_parameters.txt)

## Save the numbers of replicates and numbers of reads separately
i=$(echo $cur_line | cut -d' ' -f1)
j=$(echo $cur_line | cut -d' ' -f2)

## Run rMATS paired with the current parameter settings
rMATS_dir='Results_CEU/rMATS'
mkdir ${rMATS_dir}
mkdir ${rMATS_dir}/${i}_replicates_${j}_variance
cd ${rMATS_dir}/${i}_replicates_${j}_variance

## Run stand-alone version of rMATS
rmats_loc='rMATS-STAT'
python ${rmats_loc}/rMATS_paired.py /Simulated_Data_CEU/SimData_${i}_replicates_${j}_variance.txt ./ 4 0.05 > output_${i}_replicates_${j}_variance.txt

## Clean up house
mv rMATS_Result_P.txt rMATS_Result_P_${i}_replicates_${j}_variance.txt
