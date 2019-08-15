# PAIRADISE
## <a name="PAIRADISE"></a>What is PAIRADISE ?

PAIRADISE (**PAI**red **R**eplicate analysis of **A**llelic **DI**fferential **S**plicing **E**vents) is a method for detecting allele-specific alternative splicing (ASAS) from RNA-seq data. Unlike conventional approaches that detect ASAS events one sample at a time, PAIRADISE aggregates ASAS signals across multiple individuals in a population. By treating the two alleles of an individual as paired, and multiple individuals sharing a heterozygous SNP as replicates, PAIRADISE formulates ASAS detection as a statistical problem for identifying differential alternative splicing from RNA-seq data with paired replicates.

![PAIRADISE](https://github.com/Xinglab/PAIRADISE/blob/master/manuscript_figures/Figure1/Figure1.jpg)

## Table of Contents

- [What is PAIRADISE ?](#PAIRADISE)
- [Installation](#install)
  - [Operating system](#os)
  - [Cloning and building the PAIRADISE pipeline](#build)
  - [Components](#compon)
- [Getting started with toy example in `test_data`](#start)
- [Output folder](#input_output)
- [Contact](#contact)


## <a name="install"></a>Installation
### <a name="os"></a>Operating system
rMATS-ISO currently can only be built and run on Linux/Unix systems.

### <a name="build"></a>Cloning and building rMATS-ISO pipeline
```
git clone https://github.com/Xinglab/rMATS-ISO.git --recursive
cd rMATS-ISO
make
export PATH=$PATH:$PWD/lr2rmats/bin # To permanently modify your PATH, you need to add it to your ~/.profile or ~/.bashrc file.
```

Python 3 is required for the installnation of lr2rmats module that uses snakemaker. After the building is done, the path of `rMATS-ISO/lr2rmats/bin` needs to be added to the environment variable PATH.

### <a name="compon"></a>Components
<!---
rMATS-ISO consists the following components: lr2rmats, IsoModule, rMATS-EM, IsoClassify, and IsoPlot.
-->
rMATS-ISO consists the following components: IsoModule, rMATS-EM. More functions will be available in future versions.

IsoModule: Alternative splicing module detection.
```
python rMATS-ISO.py module --gtf --bam -o
```

rMATS-EM: Statistical test for differential splicing.
```
python rMATS-ISO.py stat --bam -o
```
They will all be automatically downloaded and built via `make` command.


## <a name="start"></a>Getting started with test example in `test_data`

<!---
python ./rMATS-ISO.py --in-gtf ./test_data/gtf/original.gtf --in-genome ./test_data/genome/genome.fa --in-long ./test_data/long_read_fa_input.list --in-short ./test_data/short_read_fa_input.list -o ./output

```
python ./rMATS-ISO.py --in-gtf ./test_data/gtf/PC3E_GS689.gtf --in-bam ./test_data/PC3E_GS689_short_read_bam_input.list -o ./output2/
```
-->

```
python rMATS-ISO.py module --gtf ./test_data/gtf/PC3E_GS689.gtf --bam ./test_data/PC3E_GS689_short_read_bam_input.list -o ./output2/
python rMATS-ISO.py stat --bam ./test_data/PC3E_GS689_short_read_bam_input.list -o ./output2/
```

## <a name="isomodule"></a>Module: Optional parameters in the splicing module detection step

Options:

         --novel-sj             allow novel splice junction in the ASM. [False]
                                   Novel splice junction here means novel combination of known splice sites.
         --novel-exon           allow novel exon in the ASM. [False]
                                   Novel exon here means novel splice sites or novel intron retention.
         --un-pair              set -u to use both proper and unproper paired mapped read. [False] (proper paired only)
         --use-multi            use both uniq- and multi-mapped reads in BAM file. [False] (uniq-mapped only)
         --anchor-len  [INT]    minimum anchor length for junction read. [1].
         --intron-len  [INT]    minimum intron length for junction read. [3]
         --exon-thres  [INT]    maximum number of exon for ASM to enumerate all possible isoforms. [10]
                                   For larger ASM, heuristic strategy will be applied. 

         --junc-cnt    [INT]    minimum allowed number of read count for known-junction. [1]
         --novel-jcnt  [INT]    minimum allowed number of read count for novel-junction. [3]
         --iso-exon-thres [INT] maximum allowed number of exons in one ASM. [50]
         --iso-cnt-thres  [INT] maximum allowed number of isoforms in one ASM. [20]

## <a name="input_output"></a>Output folder
The test output folder will contain these following sub-folders:
```
ISO_module/
EM_out/
```

The Iso_module folder contains the detected splicing modules in files ended with IsoExon. The first line of each module contains the module ID, number of exons in the module, the number of isoforms, strand, chromosome, and gene name. The second line contains the exon coordinates. The rest lines contains the isoform definations in each row, including the how the exons are included in the module and how exons are connected in the isoform. Each row is an isoform.

```
ASM#0	4	3	+	chr2	MYO1B	ENSG00000128641	A
192265108,192265194	192265475,192265561	192267358,192267444	192272841,192272915	
4	0	1	2	3
3	0	2	3
2	0	3
```

The sample-specific read counts for each module are in files ended with IsoMatrix. The first line of each module contains the module ID, the number of patterns that a read can be overlapped with multiple isoforms, the number of isoforms, the left and right side length added to the module for reads to be fully covered.

Note that for a total of N isoforms, the maximum number of overlap patterns is 2^N-1, for reads to be overlapped with different combinations of isoforms.

```
ASM#0	5	3	219	135
101	776	0	0	1
101	157	0	1	0
101	133	1	0	0
101	101	1	1	0
101	1802	1	1	1
```

The EM_out folder contains the p values for the differential splicing of modules between two sample groups and a few additional parameters from the samples. 

The p value for the differnetial splicing of splicing module is in the 6th column.

```
ASM_name total_isoforms total_exons total_read_count_group_1 total_read_count_group_2 p_value test_statistic isoform_inclusion_group_1 isoform_inclusion_group_2 isoform_inclusion_constrained variance_group_1 variance_group_2 variance_constrained dirichlet_parameter_group_1 dirichlet_parameter_group_2 dirichlet_parameter_constrained paired_isoform_pvalues isoform_index merge_info 
ASM#0 3 4 8904 6839 3.02868841117743e-13 61.347670914598 0.0686352968371457 0.575591066075565 0.277663108213901 0.000524292503279509 0.00139804133727612 0.0382402301701812 8.29974021020214 100 1.17865302356842 0.00138 1 NA 
```

## <a name="contact"></a>Contact

Levon Demirdjian demirdjial@email.chop.edu

Yi Xing xingyi@email.chop.edu

[github issues](https://github.com/Xinglab/PAIRADISE/issues)
