#!/python

# Copyright (C) 2015 University of California, Los Angeles (UCLA)
# Shayna R. Stein, Emad Bahrami-Samani, Yi Xing
#
# Authors: Shayna R. Stein, Emad Bahrami-Samani, Yi Xing
# Most recent Update: 07/19/2019 by Yungang Xu (yungag.xu@hotmail.com)
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

# imports
import sys, os, shutil,gzip
import subprocess
import re,logging,time,datetime,commands,argparse,random
from collections import defaultdict
import pysam, pybedtools
from pysam import Samfile
from pybedtools import BedTool
from  more_itertools import unique_everseen
#import time;
#ts = time.time()
## function and class definitions

#def bam2sort(prefix):
#  sys.stdout.write("sort and index: "+ prefix + '\n')
#  bam_fn = prefix + '.bam'
#  pysam.sort(bam_fn, prefix+'.sorted')
#  pysam.index(prefix+'.sorted.bam')
#  return


class AlleleAssignment :
  def __init__(self, outDir, vcf, hap1Bam, hap2Bam,  writeConflicting,rnaedit,editFile,gzipped,nrefBam,nmask,nomerge):
    self.outDir = outDir
    self.genoInfo = os.path.join(outDir,".geno") ## modified by Yungang Xu; 07/19/2019
    if os.path.isdir(vcf):
      self.vcf = vcf
      self.vcfdir = True
    else:
      self.vcf = vcf
      self.vcfdir = False
    self.hap1Bam = hap1Bam
    self.hap2Bam = hap2Bam
    self.conflicting = writeConflicting
    self.VCFheader = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']
    self.rnaedit = rnaedit
    self.editFile = editFile
    self.report = os.path.join(outDir, "report.allele_assignment.txt")
    self.gzipped = gzipped
    self.nrefBam = nrefBam
    self.nmask = nmask
    self.nomerge = nomerge

  def read_in_rna_editing(self):
    sys.stdout.write('processing: ' + self.editFile + '\n')
    e = defaultdict(lambda:defaultdict(list))
    fin = open(self.editFile)
    header = []
    firstline = True
    for line in fin:
      if firstline:
        header = line.rstrip().split()
        firstline = False
      else:
        fields = line.rstrip().split()
        result = {}
        for i,col in enumerate(header):
          result[col] = fields[i]
        e[result['chromosome'][3:]][(int(result['position'])-1)/100].append(int(result['position'])-1) 
    return e

    def read_in_vcf_nmask(self):
      if self.vcfdir:
        fList = [os.path.join(self.vcf,f) for f in os.listdir(self.vcf)  if f.endswith('.vcf') or f.endswith('.vcf.gz')]## Yungang modified: only read in the vcf or vcf.gz files
      else:
        fList = [self.vcf]
          ## output the hap info for AF calculation and sashimi plots
      self.fout_geno = open(self.genoInfo,'w') ## modified by Yungang Xu; 07/19/2019
      self.fout_geno.write('CHROM'+'\t'+'POS'+'\t'+'ID'+'\t'+'REF'+'\t'+'ALT'+'\t'+'HAP1'+'\t'+'HAP2\n')## modified by Yungang Xu; 07/19/2019
      v = defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
      vids = defaultdict(lambda:defaultdict(lambda:defaultdict(str)))
      for f in fList:
        sys.stdout.write("processing: " +f + "\n")
        if self._gzipped:
          vcf_in = gzip.open(f)
        else:
          vcf_in = open(f)
        for line in vcf_in:
          if line.startswith('#'): # skip over header lines 
            continue
          result = {} # store line in dictionary
          fields = line.rstrip().split()
          for i,col in enumerate(self.VCFheader):
            result[col] = fields[i]
          infos = [x for x in result['INFO'].split(';') if x.strip()]
          for i in infos:
            if '=' in i:
              key,value = i.split('=')
              result[key] = value
          geno = result['SAMPLE'].split(':')[0]
          alt = result['ALT'].split(',')
          alleles = [result['REF']] + alt 
          if all([re.match(r'[ACGT]',i) for i in alleles]) and all([len(i)==1 for i in alleles]):
            try: 
              g1 = int(geno.split('/')[0]) ## hap1 genotype
              g2 = int(geno.split('/')[1]) ## hap2 genotype
            except:
              g1 = int(geno.split('|')[0]) ## hap1 genotype 
              g2 = int(geno.split('|')[1]) ## hap2 genotype 
            if (g1!=g2): # heterozygous snp                                               
              v[result['CHROM']][(int(result['POS'])-1)/100][int(result['POS'])-1] = [result['REF'],result['ALT']]
              vids[result['CHROM']][(int(result['POS'])-1)/100][int(result['POS'])-1] = result['ID']
            self.fout_geno.write(result['CHROM']+'\t'+result['POS']+'\t'+result['ID']+'\t'+alleles[0]+'\t'+result['ALT']+'\t'+alleles[g1]+'\t'+alleles[g2]+'\n') ## modified by Yungang Xu; 07/19/2019
        vcf_in.close()
      self.fout_geno.close ## modified by Yungang Xu; 07/19/2019
      return v,vids

  def read_in_vcf(self):
    v = defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    vids = defaultdict(lambda:defaultdict(lambda:defaultdict(str)))
    ## output the hap info for AF calculation and sashimi plots
    self.fout_geno = open(self.genoInfo,'w') ## modified by Yungang Xu; 07/19/2019
    self.fout_geno.write('CHROM'+'\t'+'POS'+'\t'+'ID'+'\t'+'REF'+'\t'+'ALT'+'\t'+'HAP1'+'\t'+'HAP2\n')## modified by Yungang Xu; 07/19/2019
    if self.vcfdir:
      fList = [os.path.join(self.vcf,f) for f in os.listdir(self.vcf)  if f.endswith('.vcf') or f.endswith('.vcf.gz')]## Yungang modified: only read in the vcf or vcf.gz files
    else:
      fList = [self.vcf]
    n=0
    fNum=len(fList)
    for f in fList:
      n += 1
      if self.gzipped:
        vcf_in = gzip.open(f)
      else:
        vcf_in = open(f)
      sys.stdout.write("Reading in VCF:("+str(n)+"/"+str(fNum)+") " + str(f) + "\n")
      for line in vcf_in:
        if line.startswith('#'): # skip over header lines
          continue
        result = {} # store line in dictionary
        fields = line.rstrip().split()
        for i,col in enumerate(self.VCFheader):
          result[col] = fields[i]
        infos = [x for x in result['INFO'].split(';') if x.strip()]
        for i in infos:
          if '=' in i:
            key,value = i.split('=') 
            result[key] = value
        alt = result['ALT'].split(',')
        alleles = [result['REF']] + alt
        if all([re.match(r'[ACGT]',i) for i in alleles]) and all([len(i)==1 for i in alleles]):
          geno = result['SAMPLE'].split(':')[0]
          if '|' in geno: ## if genotype is phased
            g1 = int(geno.split('|')[0]) ## hap1 genotype
            g2 = int(geno.split('|')[1]) ## hap2 genotype
            if (g1!=g2): # heterozygous snp
              v[result['CHROM']][(int(result['POS'])-1)/100][int(result['POS'])-1] = [alleles[g1],alleles[g2],alleles[0]] ## subtract one from position (1 based) to match bam file (0 based)
              vids[result['CHROM']][(int(result['POS'])-1)/100][int(result['POS'])-1] = result['ID']
            else:
              vids[result['CHROM']][(int(result['POS'])-1)/100][int(result['POS'])-1] = result['ID']
            self.fout_geno.write(result['CHROM']+'\t'+result['POS']+'\t'+result['ID']+'\t'+alleles[0]+'\t'+result['ALT']+'\t'+alleles[g1]+'\t'+alleles[g2]+'\n') ## modified by Yungang Xu; 07/19/2019
      vcf_in.close()
    self.fout_geno.close ## modified by Yungang Xu; 07/19/2019
    return v,vids

  def get_tags(self,r):
    return {key: value for key, value in r.tags}
  
  def num_mismatches(self,r):
    t = {key: value for key, value in r.tags}
    return int(t['nM'])

  def haplotype_specific_read(self,r,v,h,e):
    # return 0=hap spec, 1=not hap spec, 2=conflicting, 5 = doesn't cover het snp 
    reference_positions = r.get_reference_positions()
    edit = [p for p in reference_positions if p in e[p/100]]
    snps = [p for p in reference_positions if ((p in v[p/100]) and (p not in edit))]
    if len(snps)==0: # read doesn't cover het snp
      return 5, edit, [], []
    check = [0 if r.seq[reference_positions.index(p)]==v[p/100][p][h] else 3 for p in snps]
    allele = ['R' if v[p/100][p][h]==v[p/100][p][2] else 'A' for p in snps]
    if check.count(0)>check.count(3):  # majority votes for read being hap spec
      return 0,edit,snps,allele
    elif check.count(0)==check.count(3): # conflicting
      return 2,edit,snps,allele
    else: #majority votes for other haplotype
      return 1,edit,snps,allele
      
  def haplotype_specific_read_nmask(self,r,v,e):
    # 0 = covers het snps, 3 = multimapping, 5 = no snps
   reference_positions = r.get_reference_positions()
   edit = [p for p in reference_positions if p in e[p/100]]
   snps = [p for p in reference_positions if p in v[p/100]]
   if len(snps)==0:
     return 5,edit,[],[]
   else:
     allele = ['R' if r.seq[reference_positions.index(p)]==v[p/100][p][0] else 'A' for p in snps]
     return 0,edit,snps,allele



  def haplotype_assignment_bychrom(self,chrom,hetsnps,snpids,editpositions):
    bam1 = pysam.Samfile(self.hap1Bam)
    bam2 = pysam.Samfile(self.hap2Bam)
    snpreads1,snpreads2 = defaultdict(list),defaultdict(list)
    snpreads = defaultdict(lambda:defaultdict(list))
    reads1,reads2 = defaultdict(list),defaultdict(list)

    spec1,spec2  = list(),list()
    snps1,snps2 = defaultdict(list),defaultdict(list)
    refalt1,refalt2 = defaultdict(list),defaultdict(list)
    edit1,edit2= defaultdict(list),defaultdict(list)
    hap1only,hap2only = list(),list()
    sys.stdout.write( "Reading in hap1 bam file: "+ self.hap1Bam + "\n")

    for r in bam1.fetch('chr'+str(chrom)):
      tags = self.get_tags(r)
      if int(tags['NH'])>1: ## read is multimapped, deal with separately
        continue
      spec,editpos,snppos,refalt = self.haplotype_specific_read(r,hetsnps,0, editpositions)
      snpreads1[spec].append(r.qname)
      reads1[r.qname].append(r)
      edit1[r.qname] += editpos
      snps1[r.qname] += snppos
      refalt1[r.qname] += refalt
      
    sys.stdout.write( "Reading in hap2 bam file: "+ self.hap2Bam+"\n")
    for r in bam2.fetch('chr'+str(chrom)):
      tags = self.get_tags(r)
      if int(tags['NH'])>1: ## read is multimapped
        continue
      spec,editpos,snppos,refalt = self.haplotype_specific_read(r,hetsnps,1, editpositions)
      snpreads2[spec].append(r.qname)
      reads2[r.qname].append(r)
      edit2[r.qname]+=editpos
      snps2[r.qname]+=snppos
      refalt2[r.qname] += refalt

    sys.stdout.write("Assign haplotype specific reads\n")
    conflicting = list(set([ r for r in snpreads1[0] if r in snpreads2[0]] + snpreads1[2] + snpreads2[2]))
    multimapped = list(set(snpreads1[3] + snpreads2[3]))
    rnaeditreads = list(set(snpreads1[4] + snpreads2[4]))
    snpreads1[0] = list(set([ r for r in snpreads1[0] if ((r not in conflicting) and (r not in snpreads2[3]))]))
    snpreads2[0] = list(set([ r for r in snpreads2[0] if ((r not in conflicting) and (r not in snpreads1[3]))]))
  
    for qname in snpreads1[0]:
      if qname in reads2:
        if all([True if reads1[qname][i].pos==reads2[qname][i].pos else False for i in range(len(reads1[qname]))]): # reads have same starting position in hap1 and hap2
          if self.num_mismatches(reads1[qname][0]) < self.num_mismatches(reads2[qname][0]):
            spec1.append(qname)
      else:
        hap1only.append(qname)

    for qname in snpreads2[0]:
      if qname in reads1:
        if all([True if reads1[qname][i].pos==reads2[qname][i].pos else False for i in range(len(reads2[qname]))]): # reads have same starting position in hap1 and hap2 
          if self.num_mismatches(reads2[qname][0]) < self.num_mismatches(reads1[qname][0]):
            spec2.append(qname)
      else:
        hap2only.append(qname)
    
    conflictCount = len(set(conflicting))
    hap1Count = len(set(spec1))
    hap2Count = len(set(spec2))

    out1_fn = os.path.join(self.outDir, "hap1."+str(chrom)+".bam")
    out2_fn = os.path.join(self.outDir, "hap2."+str(chrom)+".bam")
    out1 = pysam.Samfile(out1_fn,'wb',template=bam1)
    out2 = pysam.Samfile(out2_fn,'wb',template=bam2)

    sys.stdout.write("Printing bam files: " + out1_fn + ", " + out2_fn + "\n")
    for qname in spec1:
      for r in reads1[qname]:
        r.tags += [('HT',1),('SP',';'.join([str(p) for p in snps1[qname]])),('GT',';'.join([x for x in refalt1[qname]]))]
        if len(edit1[qname])>0:
          r.tags += [('EP',';'.join([str(s) for s in edit1[qname]]))]
        out1.write(r)

    for qname in spec2:
      for r in reads2[qname]:
        r.tags += [('HT',2),('SP',';'.join([str(p) for p in snps2[qname]])),('GT',';'.join([x for x in refalt2[qname]]))]
        if len(edit2[qname])>0:
          r.tags += [('EP',';'.join([str(s) for s in edit2[qname]]))]
        out2.write(r) 

    if self.conflicting:
      sys.stdout.write("Printing conflicting reads bam file \n")
      conflict1 = pysam.Samfile(self.outDir + "/hap1."+chrom+'.conflicting.bam','wb',template=bam1)
      conflict2 = pysam.Samfile(self.outDir + "/hap2."+chrom+'.conflicting.bam','wb',template=bam2)
      for qname in conflicting:
        for r in reads1[qname]:
          conflict1.write(r)
        for r in reads2[qname]:
          conflict2.write(r) 

    return hap1Count,hap2Count,conflictCount

  def haplotype_assignment(self):
    sys.stdout.write("Reading in VCF and editing files\n")
    VCF,VCFids = self.read_in_vcf()
    RNAedit = self.read_in_rna_editing()
    counts = defaultdict(lambda:defaultdict(int))
    input_files1,input_files2 = list(), list()
    for c in VCF:
      sys.stdout.write("Assigning haplotype reads in chromosome "+str(c) + '\n')
      counts[c]['h1'], counts[c]['h2'],counts[c]['c'] = self.haplotype_assignment_bychrom(c,VCF[c],VCFids[c],RNAedit[c])
      input_files1.append(os.path.join(self.outDir,"hap1."+str(c)+".bam"))
      input_files2.append(os.path.join(self.outDir,"hap2."+str(c)+".bam"))
    sys.stdout.write("Done assigning haplotype reads\nPrinting report file\n")
    if not self.nomerge:
      report_out = open(os.path.join(self.outDir,'report.assignment.txt'),'w')
      for c in counts:
        report_out.write(str(c) + '\t' + str(counts[c]['h1']) + '\t' + str(counts[c]['h2']) + '\t' + str(counts[c]['c']) + '\n')
      report_out.close()
      sys.stdout.write("Merging chromosome bam files\n")
      merge_parameters1 = ['-f',os.path.join(self.outDir,"hap1.bam")] + input_files1
      merge_parameters2 = ['-f',os.path.join(self.outDir,"hap2.bam")] + input_files2
      if len(input_files1)>1:
        pysam.merge(*merge_parameters1)
        pysam.merge(*merge_parameters2)
      else:
        os.rename(input_files1[0],os.path.join(self.outDir,"hap1.bam"))
        os.rename(input_files2[0],os.path.join(self.outDir,"hap2.bam"))
      sys.stdout.write("Sorting and indexing haplotype specific bam files - cancelled\n")
#      bam2sort(os.path.join(self.outDir,"hap1"))
#      bam2sort(os.path.join(self.outDir,"hap2"))
      sys.stdout.write("Cleaning up files - cancelled\n")
#      os.remove(os.path.join(self.outDir,"hap1.bam"))
#      os.remove(os.path.join(self.outDir,"hap2.bam"))
#      for fn in input_files1:
#        os.remove(fn)
#      for fn in input_files2:
#        os.r - cancelled 2emove(fn)
#      os.remove(os.path.join(self.outDir,"hap1.fa"))
#      os.remove(os.path.join(self.outDir,"hap2.fa"))
    sys.stdout.write("Done!\n")
    return 

  def haplotype_specific_nmask_bychrom(self,chrom,hetsnps,snpids,editpositions):
    sys.stdout.write('# assigning haplotype specific reads to N-masked genome\n')
    bam = pysam.Samfile(self._nrefBam)
    snpreads_nref = defaultdict(list)
    reads_nref = defaultdict(list)
    sys.stdout.write("# reading in VCF file\n")
    spec_nref = list()
    refalt_nref = defaultdict(list)
    snps_nref = defaultdict(list)
    edit = defaultdict(list)
    sys.stdout.write("# assigning reads\n")
    for r in bam.fetch('chr'+str(chrom)):
      tags = self.get_tags(r)
      if int(tags['NH'])>1: ## read is multimapped
        continue
      spec,editpos,snppos,refalt = self.haplotype_specific_read_nmask(r,hetsnps,editpositions)
      snpreads_nref[spec].append(r.qname)
      reads_nref[r.qname].append(r)
      refalt_nref[r.qname] += [refalt[i]  for i in range(len(refalt)) if (snppos[i] not in snps_nref[r.qname]) ]
      snps_nref[r.qname] += [snppos[i] for i in range(len(snppos)) if snppos[i] not in snps_nref[r.qname] ]
      edit[r.qname] += editpos

    snpreads_nref[0] = list(set(snpreads_nref[0]))
    sys.stdout.write("# print bam file output\n")
    out = pysam.Samfile(self._outDir + "/nmask."+chrom+'.bam','wb',template=bam)
    for qname in snpreads_nref[0]:
      for r in reads_nref[qname]:
        r.tags += [('SP',';'.join([str(p) for p in snps_nref[qname]])),('GT',';'.join([x for x in refalt_nref[qname]])), ('RF',';'.join([hetsnps[p/100][p][0] for p in snps_nref[qname]])) , ('AT',';'.join([hetsnps[p/100][p][1] for p in snps_nref[qname]]))]
        if len(edit1[qname])>0:
          r.tags += [('EP',';'.join([str(s) for s in edit[qname]]))]
        out.write(r)
    return
  
  def haplotype_assignment_nmask(self):
    sys.stdout.write("Reading in VCF and editing files\n")
    VCF,VCFids = self.read_in_vcf_nmasked()
    RNAedit = self.read_in_rna_editing()
    input_files = list()
    for c in VCF:
      sys.stdout.write("Assigning haplotype reads in chromosome "+str(c) + '\n')
      self.haplotype_assignment_nmask_bychrom(c,VCF[c],VCFids[c],RNAedit[c])
      input_files.append(os.path.join(self.outDir,"nmask."+str(c)+".bam"))
    sys.stdout.write("Done assigning haplotype reads\nPrinting report file\n")
    sys.stdout.write("Merging chromosome bam files\n")
    merge_parameters1 = ['-f',os.path.join(self.outDir,"nmask.bam")] + input_files
    if len(input_files)>1:
      pysam.merge(*merge_parameters)
    else:
      os.rename(input_files[0],os.path.join(self.outDir,"nmask.bam"))
    sys.stdout.write("Sorting and indexing haplotype specific bam files - cancelled 2\n")
#    bam2sort(os.path.join(self.outDir,"nmask"))
    sys.stdout.write("Cleaning up files - cancelled 2\n")
#    os.remove(os.path.join(self.outDir,"nmask.bam"))
#    for fn in input_files:
#      os.remove(fn)
#    os.remove(os.path.join(self.outDir,"nmask.fa"))
    sys.stdout.write("Done!\n")
    return



def main():
  parser = argparse.ArgumentParser(description='pairadise2_assign v20190720')
  parser.add_argument('-o', type=str, default="", help='output dir')
  parser.add_argument('-v', type=str, default='', help='vcf file or folder for vcf files')
  parser.add_argument('-e', type=str, default='', help='RNA editing annotation')
  parser.add_argument('--gz',action='store_true', help='the input vcf files are in the gz formate')
  parser.add_argument('--rnaedit', help='run editing analysis', action='store_true')
  parser.add_argument('--b1', type=str, default="", help='The bam file for haplotype 1')
  parser.add_argument('--b2', type=str, default="", help='The bam file for haplotype 2')
  parser.add_argument('--conflict', help='write out conflicting reads', action='store_true')
  parser.add_argument('--nmask', help='mask genome with N', action='store_true')
  parser.add_argument('--nomerge', help='merge the log out or not', action='store_true')
  args = parser.parse_args()
  print "You are running Pairadise_assign v20190720"

  outDir = args.o
  vcf = args.v
  if args.b1:
    hap1Bam = args.b1
  else:
    hap1Bam = os.path.join(outDir,"HAP1","STARalign","Aligned.sortedByCoord.out.bam")
  if args.b2:
    hap2Bam = args.b2
  else:
    hap2Bam = os.path.join(outDir,"HAP2","STARalign","Aligned.sortedByCoord.out.bam")
  nrefBam = os.path.join(outDir,"MASK","STARalign","Aligned.sortedByCoord.out.bam")
  writeConflicting = args.conflict
  rnaedit = args.rnaedit
  editFile = args.e
  gzipped = args.gz
  nmask = args.nmask

  if not args.o:
    sys.stderr.write("pairadise2: ERROR -o outdir option is required!\n\n")
    sys.exit()
  if not args.v:
    sys.stderr.write("pairadise2: ERROR -v vcf option is required!\n\n")
    sys.exit()
  nomerge = args.nomerge

  p = AlleleAssignment(outDir, vcf, hap1Bam, hap2Bam, writeConflicting,rnaedit,editFile,gzipped,nrefBam,nmask,nomerge)
  if not nmask:
    p.haplotype_assignment()
  else:
    p.haplotype_assignment_nmask()

