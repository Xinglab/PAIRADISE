# PAIRADISE <a name="PAIRADISE"></a>
Paired Replicate Analysis of Allelic Differential Splicing Events

------------

PAIRADISE (**PAI**red **R**eplicate analysis of **A**llelic **DI**fferential **S**plicing **E**vents) is a method for detecting allele-specific alternative splicing (ASAS) from RNA-seq data. Unlike conventional approaches that detect ASAS events one sample at a time, PAIRADISE aggregates ASAS signals across multiple individuals in a population. By treating the two alleles of an individual as paired, and multiple individuals sharing a heterozygous SNP as replicates, PAIRADISE formulates ASAS detection as a statistical problem for identifying differential alternative splicing from RNA-seq data with paired replicates.

<img src="https://github.com/Xinglab/PAIRADISE/blob/master/figures/Figure1/Figure1.jpg" width="1500" height="200" />

## Table of Contents

- [Installation](#install)
- [Usage](#use)
- [Output](#output)
- [Shiny App](#shiny)
- [Contact](#contact)

## <a name="install"></a>Installation

#### Dependencies ####
PAIRADISE requires the following dependencies:
- Python >= (2.7.5)
- R >= (3.2.1)
- STAR >= (2.6.0a) (can be downloaded <a href="https://github.com/alexdobin/STAR/releases/tag/2.6.0a" target="_blank">here</a>)

The folowing Python dependencies are required:
- pysam >= (0.15.2)
- pybedtools >= (0.8.0)
- more-itertools >= (5.0.0)
- numpy >= (1.16.2)

In addition, the following R packages must be installed before installing PAIRADISE:
- nloptr >= (1.2.1)
- doParallel >= (1.0.14)
- foreach >= (1.4.4)
- binom >= (1.1.1)
- ggplot2 >= (3.1.1)

#### Download and setup ####
Download PAIRADISE from github:

```git clone github.com/Xinglab/PAIRADISE/pairadise```
  
Change your working directory to "pairadise" and run the following commands to configure and install PAIRADISE:

```
./configure
make
make install
```

To test if PAIRADISE has been successfully installed, type in the command ```pairadise_personalize -h```. You should see

```
pairadise_personalize -h
usage: pairadise_personalize [-h] [--gz] [-e E] [--rnaedit] [-v V] [-o O]
                             [-r R]
                             [command [command ...]]

pairadise2

positional arguments:
  command

optional arguments:
  -h, --help  show this help message and exit
  --gz        flag denoting gzipped reads
  -e E        file containing RNA editing positions, downloaded from RADAR
  --rnaedit   flag to check for RNA editing events, must also provide an RNA
              editing file usng -e parameter
  -v V        VCF genotype directory
  -o O        output directory
  -r R        reference fasta file

```

You may need permission if want to install PAIRADISE to a root PATH. You can bypass this issue by specifying a user R library directory:

`R CMD INSTALL -l /user/R/lib src/pairadise_model/`
  
## <a name="use"></a>Usage
PAIRADISE requires the following subdirectories and input data to be in the directory where you will be performing your analysis (we'll refer to this as the 'data directory'):

1. **genome**: Contains the reference genome gtf file, fasta files, and (optionally) RNA editing information. We used the following files:

    1. hg19.fa
    2. Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
    3. Homo_sapiens.Ensembl.GRCh37.75.gtf
    4. Human_AG_all_hg19_v2.txt

2. **genotype**: Contains the genotype data in the format of vcf files, one vcf file per chromosome. The **genotype** directory 
contains one subdirectory for each sample/replicate, each containing that sample's vcf files. The subdirectory names should match the sample names. If data are biological replicates corresponding to one common sample, you only need one subdirectory.

3. **input**: Contains the fastq files for each sample/replicate. The fastq files should have the .gz extension.
4. **scripts**: Contains the scripts used to run the PAIRADISE pipeline and statistical model.

The following .txt file should also be in the data directory:

‘sample_name.txt’: Contains the sample ID’s and RNA-seq read lengths. This file contains 3 columns: the first two columns contain the sample names and the last column contains the read length. When the data are biological replicates, column 1 contains the sample IDs of the biological replicates, and column 2 contains the names of the sample.

### Running the PAIRADISE pipeline

#### Preparation:
Run the following command to create multiple .qsub files, each of which performs different stages of the PAIRADISE analysis:

```python scripts/run_preprocess.py sample_name.txt population_name path/to/data/directory/```

For example, if we are analyzing the Geuvadis CEU population and we are in the data directory, the above command becomes:

```python scripts/run_preprocess.py CEU.txt CEU ./ ```

We will continue using CEU as our running example.

------------

#### Step 1: Personalization, mapping, and assignment
Submit each of the qsub files corresponding to Step 1 of the PAIRADISE pipeline (personalization, mapping, and assignment):

``` qsub -cwd -l h_vmem=90G qsub_files/pairadise_step1_NA12843.qsub ```

The above command performs Step 1 for one sample: NA12843. You should run an analogous qsub command for each sample. 

Note: Step 1 can take a long time to run. For the CEU dataset, the typical run-time for one sample ranged from 2-5 hours.

------------

#### Step 2: Joint annotation
Once Step 1 has been completed **for every sample**, submit the qsub file corresponding to Step 2 (joint annotation):

``` qsub -cwd -l h_vmem=10G qsub_files/pairadise_step2_annotation.qsub ```


------------

#### Step 3: Counting
Once Step 2 is finished, submit each of the qsub files corresponding to Step 3 of the PAIRADISE pipeline (counting):

```qsub -cwd -l h_vmem=15G qsub_files/pairadise_step3_NA12843.qsub```

The above command performs Step 3 for one sample: NA12843. You should run an analogous qsub command for each sample.

------------

#### Step 4: Merging the counts
Once Step 3 has been completed **for every sample**, submit the qsub file corresponding to Step 4 (merging the counts):

```qsub -cwd -l h_vmem=10G qsub_files/pairadise_step4_merge.qsub```

------------

#### Step 5: Statistical modeling and visualization
Finally, we run the PAIRADISE statistical model and create plots of the significant events.

```qsub -cwd -l h_vmem=10G qsub_files/pairadise_step5_stat.qsub```

------------

## <a name="output"></a>Output

PAIRADISE outputs a table of significant ASAS events (default significance is set to FDR <= 10%). Each row of the output table corresponds to a significant ASAS event and contains the following columns:
  1. **ExonID**: Gene name, chromosome number and strand, SNP name, genomic location, reference and alternative alleles.
  2. **IJC_REF**: Exon inclusion counts for the reference allele.
  3. **SJC_REF**: Exon skipping counts for the reference allele.
  4. **IJC_ALT**: Exon inclusion counts for the alternative allele.
  5. **SJC_ALT**: Exon skipping counts for the alternative allele.
  6. **incLen**: Effective length of the exon inclusion isoform.
  7. **skpLen**: Effective length of the exon skipping isoform.
  8. **pval**: Raw (unadjusted) PAIRADISE p-value.
  9. **qval**: PAIRADISE p-value FDR adjusted using the Benjamini-Hochberg method.
  10. **IncLevel1**: Naive psi values for reference allele samples.
  11. **IncLevel2**: Naive psi values for alternative allele samples.
  12. **AvgIncLevel1**: Average psi value for reference allele samples.
  13. **AvgIncLevel2**: Average psi value for alternative allele samples.
  14. **IncLevelDifference**: Difference in average psi values between reference and alternative allele samples.
  15. **AvgTotalCount1**: Average total number of read counts for reference allele samples.
  16. **AvgTotalCount2**: Average total number of read counts for alternative allele samples.
  17. **SampleName**: Sample names.
  18. **RefAltAllele**: Reference allele and alternative allele (in the format Ref|Alt).
  19. **AF**: Allele frequencies.
 
 For visualizing the PAIRADISE results, we recommend uploading the PAIRADISE output table directly into the accompanying Shiny app. Here's an example of a significant ASAS event visualized using our Shiny app:
 
 <img src="https://github.com/Xinglab/PAIRADISE/blob/master/figures/Shiny/CEU_SCAMP3.png" width="500" height="420" />

------------

## <a name="shiny"></a>Shiny App
We have developed a data visualization tool in R Shiny for visualizing the results of PAIRADISE. The app can be accessed
<a href="http://xingshiny.research.chop.edu:3838/PAIRADISE/" target="_blank">here</a>. Have fun!


## <a name="contact"></a>Contact

Levon Demirdjian levondem@gmail.com

Yi Xing xingyi@email.chop.edu

[Report problems here](https://github.com/Xinglab/PAIRADISE/issues)

