#!/bin/bash

#$ -cwd
#$ -o Output/array_out.$JOB_ID.stdout
#$ -e Output/array_out.$JOB_ID.stderr


## This script should be called as follows:
## qsub -l h_vmem=10G -cwd -t 1-75:1 scripts/submit_R_Simulation.sh

. /etc/profile
module load R/3.5.3

## Pull the SGE_TASK_ID'th row of the parameter file as inputs
cur_line=$(sed "${SGE_TASK_ID}q;d" data/running_parameters.txt)

## Save the numbers of replicates and numbers of reads separately
i=$(echo $cur_line | cut -d' ' -f1)
j=$(echo $cur_line | cut -d' ' -f2)
k=$(echo $cur_line | cut -d' ' -f3)

## Run R job with the current parameter setting
Rscript scripts/PAIRADISE_Simulation_CEU.R $i $j $k

