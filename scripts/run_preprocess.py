##################################################
## Prepare a master script containing all of 
## the commands needed to run PAIRADISE
##
## Authors: Emad Bahrami-Samani, Levon Demirdjian
## Last Updated: 10/25/2019
## Email: demirdjial@email.chop.edu
##################################################

## Input arguments:
## 1) Text file containing the sample IDs and read lengths. This file has 3 columns:
##    Column 1: Sample names
##    Column 2: Sample names
##    Column 3: Read lengths
## 2) Name of the population to be analyzed, e.g. 'CEU'.
## 3) Working directory where the analysis will be performed

## Load necessary modules
import glob, os, sys

## Specify the working directory
pop_name = sys.argv[2]
home_dir = sys.argv[3]
input_dir = os.path.join(home_dir, "input/")

project_name = os.path.basename(sys.argv[1]).split(".")[0]
project_fn = project_name + "_samples.txt"
project_fh = open(project_fn, "w")

qsub_fn = os.path.join(home_dir, project_name + "_pairadise_pipeline.qsub")
fh = open(qsub_fn, "w")

metadata_fh = open(sys.argv[1])
all_fns = ""

## Detect extension type of fastq files
fq_ext=''
if len(glob.glob(input_dir + '*.fastq.gz')) > len(glob.glob(input_dir + '*.fq.gz')):
	fq_ext='fastq'
else:
	fq_ext='fq'

## Write pairadise_personalize, pairadise_map, pairadise_assign for each replicate separately
fh.write("###### First set of commands to process each replicate separately\n\n")
for line in metadata_fh :
  line = line.strip()
  if line == "" : continue
  fn = line.split("\t")[0]
  fid = line.split("\t")[1]
  frl = line.split("\t")[2]

  fh.write("############# " + fn + "\n")

  fh.write("pairadise_personalize -o " + project_name + "/" + fn + " -v genotype/" + fid + " -r genome/hg19.fa -e genome/Human_AG_all_hg19_v2.txt --rnaedit --gz\n")
  fh.write("pairadise_map -o " + project_name + "/" + fn + " -s input/" + fn + "_1." + fq_ext + ".gz,input/" + fn + "_2." + fq_ext + ".gz -g genome/Homo_sapiens.Ensembl.GRCh37.75.gtf -N 6 --readlength " + frl + " --gz\n")
  fh.write("samtools index " + project_name + "/" + fn + "/HAP1/STARalign/Aligned.sortedByCoord.out.bam\n")
  fh.write("samtools index " + project_name + "/" + fn + "/HAP2/STARalign/Aligned.sortedByCoord.out.bam\n")
  fh.write("pairadise_assign -o " + project_name + "/" + fn + " -v genotype/" + fid + "/ -e genome/Human_AG_all_hg19_v2.txt --rnaedit --gz\n")

  fh.write("samtools sort -o " + project_name + "/" + fn + "/hap1.sorted.bam " + project_name + "/" + fn + "/hap1.bam\n")
  fh.write("samtools index " + project_name + "/" + fn + "/hap1.sorted.bam\n")
  fh.write("samtools sort -o " + project_name + "/" + fn + "/hap2.sorted.bam " + project_name + "/" + fn + "/hap2.bam\n")
  fh.write("samtools index " + project_name + "/" + fn + "/hap2.sorted.bam\n")

  all_fns += "" + project_name + "/" + fn + "/hap1.sorted.bam," + project_name + "/" + fn + "/hap2.sorted.bam,"

  project_fh.write(project_name + "/" + fn + "\n")

  print fn

metadata_fh.close()
project_fh.close()

## Write pairadise_annotate for all replicates jointly
fh.write("###### Generate underlying annotation for all the replicates\n\n")

fh.write("mkdir ASEvents\n")
all_fns = all_fns[:-1]
fh.write("pairadise_annotate genome/Homo_sapiens.Ensembl.GRCh37.75.gtf ASEvents/fromGTF " + all_fns + " fr-unstranded " + project_name + "\n\n")

metadata_fh = open(sys.argv[1])


## Write pairadise_count for each replicate separately
fh.write("###### Second set of commands to process each replicate separately\n\n")

for line in metadata_fh :

  line = line.strip()
  if line == "" : continue
  fn = line.split("\t")[0]
  fid = line.split("\t")[1]
  frl = line.split("\t")[2]

  fh.write("############# " + fn + "\n")
  fh.write("pairadise_count -o " + project_name + "/" + fn + " --readlength " + frl + " --anchorlength 8 \n\n")

## Write pairadise_count for all replicates jointly
fh.write("###### Merging splicing counts for all the replicates\n")
fh.write("pairadise_count --merge --samples " + project_fn + "\n\n")

## Write statistical model and plotting
fh.write("###### Running statistical model\n")
fh.write("R_location=\"$(which R)\"\n")
fh.write("\"${R_location}\"script scripts/PAIRADISE_fast.R output/ASAS.SNP.SE.JunctionReadsOnly.byPair.filtered.txt results/pairadise_results.txt FALSE 0.10\n\n")

fh.write("####### Plotting significant events\n")
fh.write("\"${R_location}\"script scripts/plot_significant_events.R results/pairadise_results_filtered_FDR10.txt 1 " + pop_name + "\n\n")

fh.close()

## Make .qsub file to make job submissions easier/more organized
make_qsubs = 'sh scripts/make_qsub_files.sh ' + home_dir + ' ' + pop_name + ' ' + sys.argv[1]
os.system(make_qsubs)
