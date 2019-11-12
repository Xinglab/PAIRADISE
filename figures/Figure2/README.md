## PAIRADISE Simulations
Running the PAIRADISE simulations found in Figure 2 of the manuscript

------------

This section will outline the steps needed to run the PAIRADISE simulations found in Figure 2 of the manuscript. All of the following steps will apply to the simulation derived from the parameters of the Geuvadis CEU dataset.

The following directories and input files should be in your working directory:

1. **data**: This directory contains the list (CEU_SE_allexons_count.txt) of SE events detected on the CEU population by the PAIRADISE pipeline, as well as the list of simulation running parameters (running_parameters.txt).
2. **scripts**: This directory contains all of the necessary scripts for running the simulations.

### Step 1: Run PAIRADISE on the Geuvadis CEU data
First, we apply PAIRADISE to the Geuvadis CEU data and save the resulting model parameters for subsequent sampling:

```Rscript scripts/run_CEU_analysis.R  ```

The R workspace is saved as "Results_CEU/CEU_PAIRADISE_Results_filtered.RData" and the model parameters are saved into the text file "Results_CEU/CEU_parameters_filtered.txt". Note that only splicing events where (0.05 <= mean(psi1_hat) <= 0.95) or (0.05 <= mean(psi2_hat) <= 0.95) are retained before running PAIRADISE (events not satisfying this criterion are filtered out). 

------------

### Step 2: Generate datasets from the PAIRADISE statistical model
Next, we generate data from the PAIRADISE statistical model using the parameter estimates saved from Step 1. The datasets are saved according to the following naming scheme:

"Simulated_Data_CEU/SimData_NUMBEROFREPLICATES_replicates_VARIANCE_variance.txt"

where 

* "NUMBEROFREPLICATES" is an integer specifying the number of replicates we want in the simulation
* VARIANCE = 1,2,3 specifies if we want low (1) medium (2) or high (3) variance in the simulation. 

After generating and saving the datasets, PAIRADISE is applied to the data and the resulting workspace is saved as:

"Results_CEU/SimData_NUMBEROFREPLICATES_replicates_VARIANCE_variance.RData"

The following command performs step 2:

``` Rscript scripts/PAIRADISE_Simulation_CEU.R NUMBEROFREPLICATES VARIANCE ```

For example, to generate and analyze a dataset with 20 replicates and high variance, run the following command:

``` Rscript scripts/PAIRADISE_Simulation_CEU.R 20 3 ```

More compactly, you can submit a job array that will handle all 15 cases (5 sample sizes * 3 variance choices) at once:

``` qsub -l h_vmem=10G -cwd -t 1-15 scripts/submit_R_Simulation.sh ```

------------

### Step 3: Run the rMATS paired model
Next, we run the rMATS paired statistical model on the simulated data. We are using the stand-alone rMATS paired statistical model, found here: https://github.com/Xinglab/rMATS-STAT.

The following command will run rMATS on one simulated data file, in this case 20 replicates and low variance:

``` python rMATS_paired.py SimData_20_replicates_1_variance.txt ./ 4 0.05 > output_20_replicates_1_variance.txt ```

More compactly, you can submit a job array that will handle all 15 cases (5 sample sizes * 3 variance choices) at once (just make sure to modify the variables 'rmats_loc' and 'rMATS_dir' in scripts/submit_rMATS_paired.sh before submitting the following job):

``` qsub -l h_vmem=15G -cwd -t 1-15 scripts/run_rMATS_paired.sh ```

------------

### Step 4: Plot the AUC/TPR curves
Finally, we plot the AUC/TPR curves of PAIRADISE and the other methods. The plots are saved into "./AUC_Curves_CEU.pdf". First, we compute the AUC/TPR using the command

``` Rscript scripts/Compute_AUC_and_TPR_CEU.R simulation/working/directory```

Finally, we generate the AUC and TPR plots using the command

``` Rscript scripts/Plot_AUC_and_TPR_CEU.R directory/where/AUC/and/TPR/tables/are/located```
