# PAIRADISE
## <a name="PAIRADISE"></a>What is PAIRADISE ?

PAIRADISE (**PAI**red **R**eplicate analysis of **A**llelic **DI**fferential **S**plicing **E**vents) is a method for detecting allele-specific alternative splicing (ASAS) from RNA-seq data. Unlike conventional approaches that detect ASAS events one sample at a time, PAIRADISE aggregates ASAS signals across multiple individuals in a population. By treating the two alleles of an individual as paired, and multiple individuals sharing a heterozygous SNP as replicates, PAIRADISE formulates ASAS detection as a statistical problem for identifying differential alternative splicing from RNA-seq data with paired replicates.

![PAIRADISE](https://github.com/Xinglab/PAIRADISE/blob/master/manuscript_figures/Figure1/Figure1.jpg)

## Table of Contents

- [What is PAIRADISE ?](#PAIRADISE)
- [Installation](#install)
- [Using PAIRADISE] (#use)
- [Contact](#contact)

## <a name="install"></a>Installation
Install the dependencies:
- Python >= (2.7.5)
- R >= (3.2.1)
- STAR >= (2.6.1c)
- pysam >= (0.15.2)
- pybedtools >= (0.8.0)

In addition, the following R packages must be installed before installing PAIRADISE:
- nloptr >= (1.2.1)
- doParallel >= (1.0.14)
- foreach >= (1.4.4)

Download PAIRADISE from github:
	git clone github.com/Xinglab/PAIRADISE/pairadise
  
Change your working directory to "pairadise", `cd pairadise`, and run the following commands to configure and install PAIRADISE:

```
./configure
make
make install
```

You may need permission if want to intall PAIRADISE to a root PATH. You can bypass this issue by specifying a user R library directory:

`R CMD INSTALL -l /user/R/lib src/pairadise_model/`
  
## <a name="use"></a>Using PAIRADISE

## <a name="contact"></a>Contact

Levon Demirdjian demirdjial@email.chop.edu

Yi Xing xingyi@email.chop.edu

[github issues](https://github.com/Xinglab/PAIRADISE/issues)
