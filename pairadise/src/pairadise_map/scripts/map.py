#!/python

# Copyright (C) 2015 University of California, Los Angeles (UCLA)
# Shayna R. Stein, Emad Bahrami-Samani, Yi Xing
#
# Authors: Shayna R. Stein, Emad Bahrami-Samani, Yi Xing
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see http://www.gnu.org/licenses/.

import sys,os,gzip,re,logging,time,datetime,commands,argparse,random,shutil
from collections import defaultdict
import pysam
from pysam import Samfile

## function definitions
def STAR_create_genome(project, genome, gnme, threads,gtf,readlength):
  # build up options
  sys.stdout.write("Building genome index \n")
  opts = ""
  opts += (" --runMode genomeGenerate")
  opts += (" --genomeDir " + os.path.join(str(project),str(gnme),"STARindex"))
  opts += (" --genomeFastaFiles " + str(genome))
  opts += (" --runThreadN "+str(threads))
  opts += (" --sjdbGTFfile " + str(gtf))
  opts += (" --sjdbOverhang " + str(readlength))

  env_cpy = os.environ.copy()
  commandSTAR = ("STAR" + " " + opts)
  status,output=commands.getstatusoutput(commandSTAR)

  return

def STAR_perform_mapping(genomedir,project, gnme, seqs, threads,mismatches,gz, multimapped):
  sys.stdout.write("Aligning reads\n")
  # build up options
  opts = ""
  opts += (" --genomeDir " + str(genomedir))
  opts += (" --readFilesIn " + seqs + ' ')
  if gz:
    opts += (" --readFilesCommand gunzip -c")
  opts += (" --runThreadN " + str(threads))
  opts += (" --outFilterMultimapNmax " + str(multimapped))
  opts += (" --alignEndsType EndToEnd")
  opts += (" --outFilterMismatchNmax " + str(mismatches))
  opts += (" --outFileNamePrefix " + os.path.join(str(project), str(gnme),'STARalign/'))
  opts += (" --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated")
  opts += (" --alignIntronMax 300000 --outSJfilterOverhangMin 8 8 8 8")
  opts += (" --outSAMtype BAM SortedByCoordinate") ## output as bam
  opts += (" --outFilterMultimapScoreRange 0") ## only output top scoring reads

  env_cpy = os.environ.copy()
  commandSTAR = ("STAR" + " " + opts)
  status,output=commands.getstatusoutput(commandSTAR)

  return

## done defining functions

## begin main

def main():
  ## main function for mapping reads to personal genomes

  parser = argparse.ArgumentParser(description='pairadise2: map')
  parser.add_argument('command',nargs='*') ## command user wants to run

  ## multiple modules
  parser.add_argument('--gz',help='flag denoting gzipped reads',action='store_true') ## gzipped reads (for mapping) or VCF file (for everything else) 
  parser.add_argument('-g',help="GTF annotation file")
  parser.add_argument('-o',help="output directory")

  ## mapping parameters
  parser.add_argument('-N',help="number of mismatches per read pair")
  parser.add_argument('-T',help="num threads for STAR alignment")
  parser.add_argument('-M',help="max number of multiple alignments in STAR mapping")
  parser.add_argument('-s',help="fastq read file(s), comma deliminated if paired end")
  parser.add_argument('--readlength',help="read length, used to generate genome index (see STAR parameter sjdbOverhang)")

  args = parser.parse_args()
  command = args.command
 
  helpStr = "Help!\n"
  if not args.o:
    sys.stderr.write('pairadise_map ERROR: must provide output directory \n\n')
    sys.exit()
  if not args.s:
    sys.stderr.write('pairadise_map ERROR: must provide read sequence files \n\n')
    sys.exit()
  
  command = args.command
  if args.T:
    threads = int(args.T)
  else:
    threads = 8

  if args.N:
    mismatches = int(args.N)
  else:
    mismatches = 3

  if args.M:
    multimapped = int(args.M)
  else:
    multimapped = 20

  gzipped = args.gz
  outDir = args.o
  if not os.path.exists(outDir):
    os.makedirs(outDir)
  gtf = args.g
  hap1Ref = os.path.join(outDir, "hap1.fa")
  hap2Ref = os.path.join(outDir, "hap2.fa")

  if args.readlength:
    readlength = int(args.readlength)-1
  else:
    readlength = 99

  if len(command)>1:
    sys.stderr.write(helpStr + "\n\n")
    sys.exit()

  seqs = ' '.join((args.s).split(','))
  if len((args.s).split(','))==0 or len((args.s).split(','))>2:
    sys.stderr.write("ERROR: Sequence parameter -s input is  not correct\n Example: pairadise_map alleles -s reads_1.fq,reads_2.fq -o pairadise\n")
    sys.exit()
  
  
  if not os.path.exists(os.path.join(outDir, "HAP1/STARindex")):
    os.makedirs(os.path.join(outDir, "HAP1/STARindex"))
  if not os.path.exists(os.path.join(outDir, "HAP2/STARindex")):
    os.makedirs(os.path.join(outDir, "HAP2/STARindex"))
  if not os.path.exists(os.path.join(outDir, "HAP1/STARalign")):
    os.makedirs(os.path.join(outDir, "HAP1/STARalign"))
  if not os.path.exists(os.path.join(outDir, "HAP2/STARalign")):
    os.makedirs(os.path.join(outDir, "HAP2/STARalign"))
  if not args.g:
    sys.stderr.write('pairadise_map ERROR: must provide gtf file\n') 
    sys.exit()
  STAR_create_genome(outDir, hap1Ref, "HAP1", threads, gtf, readlength)
  STAR_create_genome(outDir, hap2Ref, "HAP2", threads, gtf, readlength)
  genomeDir1 = os.path.join(outDir, 'HAP1/STARindex')
  genomeDir2 = os.path.join(outDir, 'HAP2/STARindex')
  STAR_perform_mapping(genomeDir1, outDir, "HAP1", seqs,threads,mismatches,gzipped,multimapped)
  STAR_perform_mapping(genomeDir2, outDir, "HAP2", seqs,threads,mismatches,gzipped,multimapped)
    
