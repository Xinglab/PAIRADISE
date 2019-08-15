#
## this program processes  GTF file and SAM files to get both known and novel AS events
#

### import necessary libraries
import re,os,sys,logging,time,datetime,pickle;
import pysam
from pysam import Samfile
#
myVer= "3.2.1";
#
### checking out the number of arguments
if (len(sys.argv)<5):
  print('Not enough arguments!!');
  print ('It takes at least 4 arguments.');
  print ('Usage:\n\tpython ProgramName.py gtfFile outputPrefix SAMfiles libType [logFolder]');
  print ('Example\n\tpython ProgramName.py AceView.ncbi_37.gtf fromAceView sample1.sam,sample2.sam,... fr-unstranded temp');
  sys.exit();


def listToString(x):
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;

def uniq(inlist):
    # order preserving
    uniques = []
    for item in inlist:
        if item not in uniques:
            uniques.append(item)
    return uniques

### setting up the logging format 
logFolder = '';
if len(sys.argv)>5: ## we have a log folder
  logFolder = sys.argv[5].strip()+'/';
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=logFolder+'log.process.GTF.SAMs.'+myVer+'.'+ str(datetime.datetime.now()),
                    filemode='w')

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

###
iFile = open(sys.argv[1]); ## input gtf file
oFile_3 = open(sys.argv[2]+'.A3SS.txt', 'w'); ## alt-3 SS output file
oFile_5 = open(sys.argv[2]+'.A5SS.txt', 'w'); ## alt-5 SS output file
oFile_ce = open(sys.argv[2]+'.SE.txt', 'w'); ## skipped exon (cassette exon) output file
oFile_mxe = open(sys.argv[2]+'.MXE.txt', 'w'); ## Mutually exclusive exon output file
#oFile_afe = open(sys.argv[2]+'.AFE.txt', 'w'); ## alt first exon output file
#oFile_ale = open(sys.argv[2]+'.ALE.txt', 'w'); ## alt last exon output file
oFile_ri = open(sys.argv[2]+'.RI.txt', 'w'); ## retained intron output file
#
neFile_3 = open(sys.argv[2]+'.novelEvents.A3SS.txt', 'w'); ## alt-3 SS novel events output file
neFile_5 = open(sys.argv[2]+'.novelEvents.A5SS.txt', 'w'); ## alt-5 SS novel events output file
neFile_ce = open(sys.argv[2]+'.novelEvents.SE.txt', 'w'); ## skipped exon (cassette exon) novel events output file
neFile_mxe = open(sys.argv[2]+'.novelEvents.MXE.txt', 'w'); ## mxe novel events output file
neFile_ri = open(sys.argv[2]+'.novelEvents.RI.txt', 'w'); ## ri events output file
#neFile_afe = open(sys.argv[2]+'.novelEvents.AFE.txt', 'w'); ## afe events output file
#neFile_ale = open(sys.argv[2]+'.novelEvents.ALE.txt', 'w'); ## ale events output file
#
samFiles = sys.argv[3].split(','); ## coulbe be multiple sam files
#
libType = sys.argv[4]; ## library type

### write out header
ceHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE";
oFile_ce.write(ceHeader+'\n');
neFile_ce.write(ceHeader+'\n');

mxeHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE";
oFile_mxe.write(mxeHeader+'\n');
neFile_mxe.write(mxeHeader+'\n');
#oFile_mxe_filtered.write(mxeHeader+'\n');

altSSHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE";
oFile_3.write(altSSHeader+'\n');
oFile_5.write(altSSHeader+'\n');
neFile_3.write(altSSHeader+'\n');
neFile_5.write(altSSHeader+'\n');

#altFLHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\tdistalExonStart_0base\tdistalExonEnd\tproximalES\tproximalEE\tflankingES\tflankingEE"
#oFile_afe.write(altFLHeader+'\n');
#oFile_ale.write(altFLHeader+'\n');
#neFile_afe.write(altFLHeader+'\n');
#neFile_ale.write(altFLHeader+'\n');

riHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\triExonStart_0base\triExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE";
oFile_ri.write(riHeader+'\n');
neFile_ri.write(riHeader+'\n');

c=0;

chunk=1000;
geneGroup={}; ## genes in a group
jGroup={}; ## junctions in a group

genes = {}; 
supple = {}; 
junctions={};  ## junctions by sample
tJunctions={}; ## total junctions

mil=50;  ## minimum intron length for novel exon or splice site finding
mel=250; ## maximum exon length for novel exon or splice site finding


##cds={}; ## coding region

for line in iFile: ## for each line
  if line.strip()[0]=='#': ## comments, skip this line
    continue;
  ele = line.strip().split('\t');
  chr = ele[0];
  type = ele[2].lower(); ## exon, intron, CDS, start_codon, stop_codon..
  sC = ele[3]; ## start coord, 1-base
  eC = ele[4]; ## end coord, 1-base
  group = range(int(sC)/chunk, int(eC)/chunk + 1); ## groups this line could belong to
  group = list(set(group));  ## remove duplicate groups
  strand = ele[6]; 
  desc = ele[8].split(';');
  gID=['','']; txID=['','']; gName=['','NA'] ## init..

  for dEle in desc: ## for each element of description
    if len(dEle.strip())<2 or len(dEle.strip().split(' '))<2:
      continue; ## probably the last description
    dName = dEle.strip().split(' ')[0];
    dVal = dEle.strip().split(' ')[1];
    if dName.upper() == 'GENE_ID': ## it is a description for gene_id
      gID = [dName,dVal];
    elif dName.upper() == 'TRANSCRIPT_ID': ## it is a description for transcript_id
      txID = [dName, dVal];
    elif dName.upper() == 'GENE_NAME': ## it is a description for gene_name
      gName = [dName, dVal];

  if gID[0].upper()!='GENE_ID' or txID[0].upper() != 'TRANSCRIPT_ID': ## wrong one..
    logging.debug("This line does not have correct description for gID or txID: %s, %s" % (gID, txID));
    logging.debug("Incorrect description: %s" % ele);
    continue; ## process next line

  for i in group: ## for each possible group
    if i in geneGroup: ## this group already exist
      geneGroup[i].append(gID[1]); ## duplicate geneIDs will get removed after the outer for loop
    else: ## first time accesing this group
      geneGroup[i] = [gID[1]];

  if type=='exon':  ## process exon
    if gID[1] in genes: # already processed this gID
      if txID[1] in genes[gID[1]]: ## this transcript is added already
        genes[gID[1]][txID[1]].append([int(sC), int(eC)]); ## add exon to the existing Tx
      else: ## first time processing this Tx
        genes[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first exon
    else:  ## new gene ID
      genes[gID[1]] = {};
      genes[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first exon
      #supple[gID[1]] = [gID[1], chr, strand]; ## geneName, chromosom and strand
      supple[gID[1]] = [gName[1], chr, strand]; ## geneName, chromosom and strand
#  if type=='cds': ## coding region
#    if gID[1] in cds: # already processed this gID
#      if txID[1] in cds[gID[1]]: ## this transcript is added already
#        cds[gID[1]][txID[1]].append([int(sC), int(eC)]); ## add CDS to the existing Tx
#      else: ## first time processing this Tx
#        cds[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first CDS
#    else:  ## new gene ID
#      cds[gID[1]] = {};
#      cds[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first exon
#logging.debug("Done populating genes and cds dictionaries");
#
#print 
#sys.exit(0);
for gg in geneGroup: ## for all groups in geneGroup
  geneGroup[gg] = list(set(geneGroup[gg]));
#
## stats
#
logging.debug("======== stats from genes =========");
#
nGene=len(genes); ## number of genes in genes dict
nTx=0; ## number of transcripts
oneTx=0; ## number of one-tx genes
nExon = 0; ## number of exons
oneExon=0; ## number of one-exon transcripts
#
oneTxOneExon=0; ## number of one-tx genes with only one exon
#
for id in genes: ## for each gene
  nTx += len(genes[id]); 
  if len(genes[id])==1:
    oneTx += 1; ## one-transcript gene
  for tx in genes[id]: ## for each tx
    nExon += len(genes[id][tx]);
    if len(genes[id][tx])==1: ## one exon tx
      oneExon += 1;
      if len(genes[id])==1: ## one tx gene
        oneTxOneExon+=1;

logging.debug("There are %d distinct gene ID in the gtf file" % nGene);
logging.debug("There are %d distinct transcript ID in the gtf file" % nTx);
logging.debug("There are %d one-transcript genes in the gtf file" % oneTx);
logging.debug("There are %d exons in the gtf file" % nExon);
logging.debug("There are %d one-exon transcripts in the gtf file" % oneExon);
logging.debug("There are %d one-transcript genes with only one exon in the transcript" % oneTxOneExon);
if (nGene>0): ## to avoid divided by zero exception
  logging.debug("Average number of transcripts per gene is %f" % (float(nTx)/nGene));
if (nTx>0): ## to avoid divided by zero exception
  logging.debug("Average number of exons per transcript is %f" % (float(nExon)/nTx));
if (nTx-oneExon)>0: ## to avoid divided by zero exception
  logging.debug("Average number of exons per transcript excluding one-exon tx is %f" % (float(nExon-oneExon)/(nTx-oneExon)));
#
#logging.debug("======== stats from cds =========");
#
#nGene=len(cds); ## number of genes in cds dict
#nTx=0; ## number of transcripts
#oneTx=0; ## number of one-tx genes
#nExon = 0; ## number of exons
#oneExon=0; ## number of one-exon transcripts
#
#oneTxOneExon=0; ## number of one-tx genes with only one exon
#
#for id in cds: ## for each gene
#  nTx += len(cds[id]);
#  if len(cds[id])==1:
#    oneTx += 1; ## one-transcript gene
#  for tx in cds[id]: ## for each tx
#    nExon += len(cds[id][tx]);
#    if len(cds[id][tx])==1: ## one exon tx
#      oneExon += 1;
#      if len(cds[id])==1: ## one tx gene
#        oneTxOneExon+=1;
#
#logging.debug("There are %d distinct gene ID in the gtf file" % nGene);
#logging.debug("There are %d distinct transcript ID in the gtf file" % nTx);
#logging.debug("There are %d one-transcript genes in the gtf file" % oneTx);
#logging.debug("There are %d exons in the gtf file" % nExon);
#logging.debug("There are %d one-exon transcripts in the gtf file" % oneExon);
#logging.debug("There are %d one-transcript genes with only one cds in the transcript" % oneTxOneExon);
#if (nGene>0): ## to avoid divided by zero exception
#  logging.debug("Average number of transcripts per gene is %f" % (float(nTx)/nGene));
#if (nTx>0): ## to avoid divided by zero exception
#  logging.debug("Average number of cds per transcript is %f" % (float(nExon)/nTx));
#if (nTx-oneExon)>0: ## to avoid divided by zero exception
#  logging.debug("Average number of cds per transcript excluding one-cds tx is %f" % (float(nExon-oneExon)/(nTx-oneExon)));
#
logging.debug("======== stats from geneGroup =========");
#
tgi = 0;## total geneIDs
## 
for gg in geneGroup: 
  tgi += len(geneGroup[gg]);
logging.debug("There are total of %d groups and %d genes in geneGroup" % (len(geneGroup), tgi));
logging.debug("The average number of genes in each group is %f" % (float(tgi)/len(geneGroup))); 
#
logging.debug("==========================================\n");
#
#
### sort transcripts ###
#
for gID in genes:
  for tID in genes[gID]: ## sort each transcript
    if len(genes[gID][tID])==1: ## only one exon, skip it
      continue; ## next transcript..
    genes[gID][tID] = sorted(genes[gID][tID]); ## sort each transcript
#
#
#### now process SAM files to populate junction dictionaries ####
#
samIndex=0;
for s1 in samFiles: ## for each samFiles
  if len(s1.strip())<1: ## incorrect split. a user might accidently put comma at the end of input sam file list
    continue; ### just skip this entry, probably the last one though
  logging.debug("processing %s" % s1.strip());
  sFile = pysam.Samfile(s1.strip(),'rb')
  for read in sFile.fetch():
    chr = sFile.getrname(read.reference_id)
    mc = read.pos+1
    mString = read.cigarstring
#  sFile = open(s1.strip()); ## open sam file
#  for line in sFile: ## process each line
#    if len(line.strip().split('\t'))<5 or line.strip()[0]=='#' or line.strip()[0]=='@' : ## blank line or comment
#      continue;  ## go to next line
#    ele = line.strip().split('\t');
#    chr = ele[2];
    if chr[0:3]!='chr':
      chr = 'chr'+chr;
#    mc = 0;
#    try:
#      mc=int(ele[3]); ## 1 base, mapping coordinate
#    except:
#      logging.debug("Could not convert 4th element to integer in %s" % line.strip());
#      continue;
#    mString = ele[5]; ## mapping string, 50M or aMbNcMdNeMfNgM format, CIGAR string
    group = mc/chunk; ## group does not change, it's okay to check only one group for a junction read
    if 'D' in mString or 'I' in mString or 'S' in mString or 'H' in mString or 'P' in mString or 'X' in mString or '=' in mString: ## skip
      continue; ## go to next line

    ### check to see if the line is either exonic read or junction read
    split_mString = mString.split('M');
    tor = 0; ## type of read, 0 nothing, 1 exonic read, 2 junction read
    if len(split_mString)==2:
      #tor = 1; ## exonic read ##
      continue; ## go to next line
    elif len(split_mString)>=3:    ###### junction read ###########
      #tor = 2; ## junction read
      jS=mc; jE=mc-1; 
      for ec in range(0,len(split_mString)-2): ## for each coordinate
        secondNumber = int(split_mString[ec].split('N')[-1]); 
        jS = jE+secondNumber; ## 1-base
        jE = jS+int(split_mString[ec+1].split('N')[0]); ## 0-base
        key = chr+'_'+str(jS)+'_'+str(jE);

        jGrp=range(jS/chunk, jE/chunk+1);
        jGrp=list(set(jGrp));  ## remove duplicate group
        for grp in jGrp: ## for each possible group
          if grp in tJunctions:
            if key in tJunctions[grp]:
              tJunctions[grp][key] += 1;
            else:
              tJunctions[grp][key]=1;
          else:
            tJunctions[grp]={};
            tJunctions[grp][key]=1;

          if grp in junctions:
            if samIndex in junctions[grp]:
              if key in junctions[grp][samIndex]:
                junctions[grp][samIndex][key]+=1;
              else:
                junctions[grp][samIndex][key]=1;
            else: 
              junctions[grp][samIndex]={};
              junctions[grp][samIndex][key]=1;
          else:
            junctions[grp]={};
            junctions[grp][samIndex]={};
            junctions[grp][samIndex][key]=1;

  samIndex+=1; ## for the junctions dictionary index

logging.debug("Done populating junctions and tJunctions dictionary");
logging.debug("Dumping two dictionaries");
tjFile=open(logFolder+'totalJunctions.pickle','w');
jpsFile=open(logFolder+'junctions.per.sample.pickle','w');
pickle.dump(tJunctions, tjFile);
pickle.dump(junctions, jpsFile);
logging.debug("Done dumping two dictionaries");


#
### now process each gene and tJunction to add novel transcripts constructed from novel junctions, novel splice site, or novel exons 
### novel junction connects the existing exons (within a gene) but the junction is not defined in the gtf
### novel splice sites are splice sites not annotated in the gtf
### novel exons are exons not annotated in the gtf
#

def getGeneGroups(gene_id):
  tmin=1000000000; tmax=-1;
  for gtx in genes[gene_id]: ## for all tx
    if genes[gene_id][gtx][0][0]<tmin: ## new min
      tmin=genes[gene_id][gtx][0][0];
    if genes[gene_id][gtx][-1][1]>tmax: ## new max
      tmax=genes[gene_id][gtx][-1][1];
  rVal=range(tmin/chunk, tmax/chunk +1);
  return rVal;

def getJunctionType(jjj, cgg):
### now get the type of junction
### Junction type
## -1: not belong to this gene
### 0: known junction
### 1: novel junction
### 2: left anchored novel splice site
### 3: right anchored novel splice site

  jEle=jjj.split('_');
  j_chr=jEle[0];

  for i in range(1,len(jEle)-2): ## if chr contains '_'
    j_chr += '_'+jEle[i]

  if j_chr != supple[cgg][1]: ## different chromosome. return -1;
    return -1;  

  #jS = int(jEle[1]); jE = int(jEle[2]);
  jS = int(jEle[-2]); jE = int(jEle[-1]);  ## in case chromosome contains '_' in it

  tempType=100; ## temp type
  for txx in genes[cgg]: ## for each transcript, look for the smallest type, when found type 0, just return it
    cexons = genes[cgg][txx]; ## candidate exons
    uInd=-1; dInd=-1; ## upstream exon index and downstream exon index

    for ci in range(0, len(cexons)): ## examine each candidate exon
      if cexons[ci][1]==jS: ## this exon has the same end
        uInd=ci;
      if cexons[ci][0]==(jE+1): ## this exon has the same start
        dInd=ci;

    if uInd>-1 and dInd>-1: ## we have both exons here.., type 0 or type1
      if dInd-uInd == 1: ## known junction.. type 0
        return 0; 
      elif dInd-uInd>1: ## novel junction.. 
        if tempType>1: ## other temp type
          tempType=1;
        else: ## do nothing, keep going
          pass;

    elif uInd>-1 and dInd==-1: ## left anchored, type 2
      if tempType<2: ## we have better type
        pass;
      elif tempType==3: ## right anchor was already found, set it to type 1
        tempType=1;
      elif tempType>2: ## set it to tempType 2
        tempType=2; 

    elif uInd==-1 and dInd>-1: ## right anchored, type 3
      if tempType<2: ## we have better type
        pass;
      elif tempType==2: ## left anchor was already found, set it to type 1
        tempType=1;
      elif tempType>3: ## set it to tempType 2
        tempType=3;
  if tempType<100:
    return tempType;
  else: ## no type was found. in valid junction
     return -1;


numJ =0;numELA=0;numERA=0;numEB=0;

def makeNovelTX(jtt,cgg): ## make novel tx using special junction types
  global numJ,numELA,numERA,numEB; ## global counter for number of new transcripts from each type
  rValue={};
  J=1; E_LA=2; E_RA=3; E_B=4; ## novel junction (1), novel exon from left anchor(2), right anchor(3), from both (4)
  rValue[J]={}; rValue[E_LA]={}; rValue[E_RA]={}; rValue[E_B]={};

  if 1 in jtt: ### there are novel junctions
    for jcoord in jtt[1]:
      #jele = jcoord.strip().split('_'); jS = int(jele[1]); jE=int(jele[2]); 
      jele = jcoord.strip().split('_'); jS = int(jele[-2]); jE=int(jele[-1]); 

      for txx in genes[cgg]: ## for each transcript, look for the exons with at least one exons between the junction
        cexons = genes[cgg][txx]; ## candidate exons
        uInd=-1; dInd=-1; ## upstream exon index and downstream exon index
    
        for ci in range(0, len(cexons)): ## examine each candidate exon
          if cexons[ci][1]==jS: ## this exon has the same end
            uInd=ci;
          if cexons[ci][0]==(jE+1): ## this exon has the same start
            dInd=ci;
    
        if uInd>-1 and dInd>-1: ## we have both exons here.., type 0 or type 1
          if dInd-uInd>1: ## novel junction..
            key=str(cexons[uInd][0])+':'+str(cexons[uInd][1])+':'+str(cexons[dInd][0])+':'+str(cexons[dInd][1]); ## novel junction
            nex = [cexons[uInd],cexons[dInd]]; ## novel exons
            if uInd>0: ## not the first exon
              key = str(cexons[uInd-1][0])+':'+str(cexons[uInd-1][1])+':'+key; ## prev exon + novel junction
              nex = [cexons[uInd-1]]+nex;
            if dInd<len(cexons)-1: ## not the last exon
              key = key + ':'+str(cexons[dInd+1][0])+':'+str(cexons[dInd+1][1]); ## novel junction + next exon
              nex = nex+[cexons[dInd+1]];
            rValue[J][key] = nex;  ## it's okay to overwrite
          else: ### should not be here..
            logging.debug("Check this out..gene: %s, tx: %s, uExon: %s, dExon: %s" % (cg,txx,cexons[uInd], cexons[dInd]));
      ## end of for txx in genes[cg]
      ## logging.debug("Number of novel junctions: %d" % len(rValue[J]);
      for novelT in rValue[J]: ## for all novel transcript from junction
        txName = cg+'.novelJunction_'+novelT;
        numJ += 1;
        genes[cg][txName] = rValue[J][novelT];


  if 2 in jtt: ## novel splice site with left anchor, right side is unknown (could be E_LA or E_B)
    ### find the closest splice 5 prime site it's either known (E_LA) or novel (E_B)   
    
    for jcoord in jtt[2]:
      jele=jcoord.strip().split('_'); jS = int(jele[-2]); jE=int(jele[-1]); ## jS is known, jE is unknown
      if jE-jS<mil: ## intron size does not meet the minimum intron length
        continue; ## do nothing 
      exonStart=-1; ## need to get exon start coord

      for txx in genes[cgg]: ## lets find the exon start coord
        cexons=genes[cgg][txx]; ## candidate exons
        for ci in range(0,len(cexons)): ## examine each candidate exon
          if cexons[ci][1]==jS: ## it is a junction from this exon
            exonStart=cexons[ci][0]; ## found it, let's break;
            break;
        if exonStart>0: ## found it
          break;
      if exonStart==-1: ### something went wrong. log it then proceed
        logging.debug("No exon start is available for %s" % jcoord);    
        continue; ## go to the next junction  

      ## check known splice site first (E_LA)
      for txx in genes[cgg]: ## for each transcript, look for the exons whose eE is greater than jE
        cexons=genes[cgg][txx]; ## candidate exons
        for ci in range(0,len(cexons)): ## examine each candidate exon
          eLength=cexons[ci][1]-jE;          
          if eLength>0 and eLength<mel: ## it meets max exon size limit, add a new transcript
            key=str(exonStart)+':'+str(jS)+':'+str(jE+1)+':'+str(cexons[ci][1]);
            nex=[[exonStart,jS],[jE+1,cexons[ci][1]]]; ## two exons first
            if ci<len(cexons)-1: ## not the last exon
              key = key + ':'+str(cexons[ci+1][0])+':'+str(cexons[ci+1][1]); ## novel exon with left anchor
              nex = nex+[cexons[ci+1]];
            rValue[E_LA][key]=nex; ## it's okay to overwrite
      ### end of for txx
      for novelT in rValue[E_LA]: ## for all novel exons from left anchored junctions
        txName = cg+'.novelLeftAnchorExon_'+novelT;
        numELA += 1;
        genes[cg][txName] = rValue[E_LA][novelT];

      ## check novel splice site (E_B)
      if 3 not in jtt: ### do not do this part. just continue
        continue;
      for rj in jtt[3]: 
        rjel = rj.strip().split('_'); rjS=int(rjel[-2]); rjE=int(rjel[-1]); ## rjS is unknown, rjE is known
        if rjE-rjS<mil: ## intron size is too small
          continue; ## next junction please

        eLength = rjS-jE;
        if eLength>0 and eLength<mel: ## it meets max exon size limit, add a new transcript

          exonEnd=-1; ## need to get exon end coord
  
          for txx in genes[cgg]: ## lets find the exon end coord
            cexons=genes[cgg][txx]; ## candidate exons
            for ci in range(0,len(cexons)): ## examine each candidate exon
              if cexons[ci][0]==rjE+1: ## it is a junction from this exon
                exonEnd=cexons[ci][1]; ## found it, let's break;
                break;
            if exonEnd>0: ## found it
              break;
          if exonEnd==-1: ### something went wrong. log it then proceed
            logging.debug("No exon end is available for %s" % rj);
            continue; ## go to the next junction

          key=str(exonStart)+':'+str(jS)+':'+str(jE+1)+':'+str(rjS)+':'+str(rjE+1)+':'+str(exonEnd);
          nex=[[exonStart,jS],[jE+1,rjS], [rjE+1,exonEnd]]; ## three exons
          rValue[E_B][key]=nex; ## it's okay to overwrite
      ### end of for txx
#      for novelT in rValue[E_B]: ## for all novel exons from both left and right anchored junctions
#        txName = cg+'.novelExonBoth__'+novelT;
#        numEB += 1;
#        genes[cg][txName] = rValue[E_B][novelT];
    
  if 3 in jtt: ## novel splice site with right anchor, left side is unknown (could be E_RA or E_B)
    ### find the closest splice 3 prime site it's either known (E_LA) or novel (E_B)

    for jcoord in jtt[3]:
      jele=jcoord.strip().split('_'); jS = int(jele[-2]); jE=int(jele[-1]); ## jS is unknown, jE is known
      if jE-jS<mil: ## intron size does not meet the minimum intron length
        continue; ## do nothing
      exonEnd=-1; ## need to get exon end coord

      for txx in genes[cgg]: ## lets find the exon end coord
        cexons=genes[cgg][txx]; ## candidate exons
        for ci in range(0,len(cexons)): ## examine each candidate exon
          if cexons[ci][0]==jE+1: ## it is a junction from this exon
            exonEnd=cexons[ci][1]; ## found it, let's break;
            break;
        if exonEnd>0: ## found it
          break;
      if exonEnd==-1: ### something went wrong. log it then proceed
        logging.debug("No exon end is available for %s" % jcoord);
        continue; ## go to the next junction

      ## check known splice site first (E_RA)
      for txx in genes[cgg]: ## for each transcript, look for the exons whose eS is smaller than jS
        cexons=genes[cgg][txx]; ## candidate exons
        for ci in range(0,len(cexons)): ## examine each candidate exon
          eLength=jS-cexons[ci][0]+1;
          if eLength>0 and eLength<mel: ## it meets max exon size limit, add a new transcript
            key=str(cexons[ci][0])+':'+str(jS)+':'+str(jE+1)+':'+str(exonEnd);
            nex=[[cexons[ci][0],jS],[jE+1,exonEnd]]; ## two exons first
            if ci>0: ## not the first exon
              key = str(cexons[ci-1][0])+':'+str(cexons[ci-1][1])+':'+key; ## novel exon with right anchor
              nex = [cexons[ci-1]]+nex;
            rValue[E_RA][key]=nex; ## it's okay to overwrite
      ### end of for txx
      for novelT in rValue[E_RA]: ## for all novel exons from right anchored junctions
        txName = cg+'.novelRightAnchorExon_'+novelT;
        numERA += 1;
        genes[cg][txName] = rValue[E_RA][novelT];

      ## check novel splice site (E_B)
      if 2 not in jtt: ### skip this
        continue;
      for lj in jtt[2]:
        ljel = lj.strip().split('_'); ljS=int(ljel[-2]); ljE=int(ljel[-1]); ## ljS is known, ljE is unknown
        if ljE-ljS<mil: ## intron size is too small
          continue; ## next junction please

        eLength = jS-ljE;
        if eLength>0 and eLength<mel: ## it meets max exon size limit, add a new transcript

          exonStart=-1; ## need to get exon start coord

          for txx in genes[cgg]: ## lets find the exon start coord
            cexons=genes[cgg][txx]; ## candidate exons
            for ci in range(0,len(cexons)): ## examine each candidate exon
              if cexons[ci][1]==ljS: ## it is a junction from this exon
                exonStart=cexons[ci][0]; ## found it, let's break;
                break;
            if exonStart>0: ## found it
              break;
          if exonStart==-1: ### something went wrong. log it then proceed
            logging.debug("No exon start is available for %s" % lj);
            continue; ## go to the next junction

          key=str(exonStart)+':'+str(ljS)+':'+str(ljE+1)+':'+str(jS)+':'+str(jE+1)+':'+str(exonEnd);
          nex=[[exonStart,ljS],[ljE+1,jS], [jE+1,exonEnd]]; ## three exons
          rValue[E_B][key]=nex; ## it's okay to overwrite
      ### end of for txx
  ## take care of E_B here, it can be done from both junction type 2 and 3
  for novelT in rValue[E_B]: ## for all novel exons from both left and right anchored junctions
    txName = cg+'.novelExonBoth__'+novelT;
    numEB += 1;
    genes[cg][txName] = rValue[E_B][novelT];

  ##return rValue;



for cg in genes: ## for each gene ID
  egrp = getGeneGroups(cg); ## gets the groups to examine
  ##print "hi",egrp; sys.exit(-10);
  numIter=3;
  for nn in range(numIter):  ## do the iteration 
    tempJ={}; ## to avoid examining the same junction multiple times from diferent group
    jType={}; ## store junctions by type
    for eg in egrp: ## now examine junction
      if eg not in tJunctions: ## no junction in this group
        continue;
      ##print len(tJunctions[eg]); sys.exit(-10);
      for jj in tJunctions[eg]: ## for each junction in each group
        if jj in tempJ: ## already processed this junction. skip it
          continue;
        tempJ[jj]=1; ## indicate that we processed this key chr_jS_jE
        jt = getJunctionType(jj,cg); ## returns junction type
        if jt>0: ## special junction
          if jt in jType: 
            jType[jt].append(jj);
          else:
            jType[jt]=[jj];
    ## now jType has all junctions we need to process to make novel transcripts from novel junction, novel ss, and novel exon
    makeNovelTX(jType,cg); ## detect and add new transcripts 
    
        ## end of for cg in geneGroup[group]
      ## end of for each coordinate
    ## end of elif junction read
#
#
def fullyContainedInInternalExon(myExon, myGeneID): ## return true if any transcript has an internal exon fully containing the exon
  for myTxID in genes[myGeneID]:
    for intExon in genes[myGeneID][myTxID][1:-1]: ## for each internal exon
      if intExon[0]<=myExon[0] and intExon[1]>=myExon[1]: ## fully contained
        return True;
  return False;
#
logging.debug("Process each gene from dictionary");
#
numSkippingEvents=0;
numMXEvents=0;
numSSEvents=0;
num3=0;
num5=0;
numAFE=0;
numALE=0;
numRI=0; 

sEvents={};
mxEvents={};
ss3Events={};
ss5Events={};
mxEvents_filtered={}; 
afeEvents={}; ## alternative first exon
aleEvents={}; ## alternative last exon
riEvents={};  ## retained intron 

dupSE=0;
dupMXE=0;
dupSS3=0;
dupSS5=0;
dupRI=0;

filteredMXE=0;
#
for gID in genes:  ## process each gene
  supInfo = supple[gID]; ## supplementary info
  if len(genes[gID])==1: ## only one transcript, alt SS event is imposible
    continue; ## go to the next geneID
  else: ## look for alt SS event
    de={}; ## distinct exons
    for tID in genes[gID]: ## sort each transcript, merge and get distinct exons
      if len(genes[gID][tID])==1: ## only one exon, skip it 
        continue; ## next transcript..
      genes[gID][tID] = sorted(genes[gID][tID]); ## sort each transcript
      for exon in genes[gID][tID]:
        de[exon[0], exon[1]] = 1; ## it's okay to overwrite

#    dc={}; ## distinct cds
#    if gID in cds: ## this gene has cds 
#      for tID in cds[gID]: ## sort each transcript, merge and get distinct cds
#        if len(cds[gID][tID])==1: ## only one cds, skip it
#          continue; ## next transcript..
#        cds[gID][tID] = sorted(cds[gID][tID]); ## sort each transcript
#        for exon in cds[gID][tID]: ## for each cds
#          dc[exon[0], exon[1]] = 1; ## it's okay to overwrite

    ## now we have sorted transcripts and distinct exons
    ## examine each exon in de to see if the exon is involved in any types of AS events
  
    for ce in de: ## for each exon in distinct exon dictionary
      uf=[]; ## upstream flanking exons
      df=[]; ## downstream flanking exons
      for tID in genes[gID]: ## examine each transcript to see if it contains the given exon
        if [ce[0],ce[1]] in genes[gID][tID]: ## this exon is in the transcript     
          eInd = genes[gID][tID].index([ce[0],ce[1]]);
          if (0<eInd): ## it definitely has upstream flanking exon
            uf.append(genes[gID][tID][eInd-1]);
          if (eInd<len(genes[gID][tID])-1): ## it definitely has downstream flanking exon
            df.append(genes[gID][tID][eInd+1]);

      ## getting uniq upstream flanking exons and downstream flanking exons
      uf = uniq(uf);
      df = uniq(df);

      #### exon skipping (cassette exon) events ###
      for i in range(0,len(uf)): ### going through the upstream flanking exons
        f1 = uf[i]; ## first flanking exon
        for j in range(0,len(df)): ## going through the downstream falnking exons to see if there is skipping junction
          f2 = df[j];
          for tID in genes[gID]: ## examine each transcript, to see if it has an edge connecting f1 and f2
            if len(genes[gID][tID])<2: ## less than two exons, skip it
              continue;

########################## not requiring exact falnking exons ###########
            for i in range(0,len(genes[gID][tID])-1): ## for each exon in the tx
              e_1 = genes[gID][tID][i];
              e_2 = genes[gID][tID][i+1]; 
              if e_1[1]==f1[1] and e_2[0]==f2[0]: ## this tx has an edge connecting f1 and f2 but does not have ce
                key=supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(f1[1])+':'+str(f2[0]-1); ## target exon and skipping junction
                if key in sEvents: ## already have this skipping events
                  dupSE +=1;
                  continue; ## next transcript
                else: ## new key, write it
                  sEvents[key]=1;
                  numSkippingEvents += 1;
                  oFile_ce.write(str(numSkippingEvents)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(f2[0]-1)+'\t'+str(f2[1])+'\n');
                  if (tID.find("novel")>=0): #It involves novel junction
                    neFile_ce.write(str(numSkippingEvents)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(f2[0]-1)+'\t'+str(f2[1])+'\n');
###################  

####   requiring exact flanking exons ##################
#            if f1 in genes[gID][tID] and f2 in genes[gID][tID]: ## this transcript has f1 and f2
#              if genes[gID][tID].index(f1)+1 == genes[gID][tID].index(f2): ### there is an edge connecting f1 and f2
#                #key=supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(f1[0]-1)+':'+str(f1[1])+':'+str(f2[0]-1)+':'+str(f2[1]);
#                key=supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(f1[1])+':'+str(f2[0]-1); ## target exon and skipping junction
#                if key in sEvents: ## already have this skipping events
#                  dupSE +=1;
#                  continue; ## next transcript
#                else: ## new key, write it
#                  sEvents[key]=1;
#                  numSkippingEvents += 1; 
#                  oFile_ce.write(str(numSkippingEvents)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(f2[0]-1)+'\t'+str(f2[1])+'\n');
##########################################################################

      ##### mutually exclusive events #####
      for i in range(0,len(uf)): ### going through the upstream flanking exons
        f1 = uf[i]; ## first flanking exon
        for j in range(0,len(df)): ## going through the downstream falnking exons to see if there is skipping junction
          f2 = df[j];
          for tID in genes[gID]: ## examine each transcript
            if [ce[0], ce[1]] in genes[gID][tID] or f1 not in genes[gID][tID] or f2 not in genes[gID][tID]: ## ce in, f1 or f2 is not in..do not examine
              continue; ### go to next transcript in genes[gID]
            else: ## this transcript does not have ce and has both f1 and f2, let's take a look
              fromF1 = genes[gID][tID].index(f1);
              fromF2 = genes[gID][tID].index(f2);
              mxe = genes[gID][tID][fromF1+1]; ## candidate mxe
              if (fromF1+1 == fromF2-1) and (mxe[0]>ce[1]): ### this exon is the right one
                goodMXE=True;
#                for txs in genes[gID]: ### search through transcripts again to see if there is ce-mxe edge
#                  if [ce[0],ce[1]] in genes[gID][txs] and mxe in genes[gID][txs]: ## ce and mxe are not MXE
#                    goodMXE = False;
#                    filteredMXE+=1;
#                    fkey = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(mxe[0]-1)+':'+str(mxe[1])+':'+str(f1[0]-1)+':'+str(f1[1])+':'+str(f2[0]-1)+':'+str(f2[1]);
#                    if fkey in mxEvents_filtered: ## duplicate, do not write out
#                      pass;
#                    else:
#                      mxEvents_filtered[fkey] = 1;
#                      oFile_mxe_filtered.write(str(filteredMXE)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(mxe[0]-1)+'\t'+str(mxe[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(f2[0]-1)+'\t'+str(f2[1])+'\n');
#                    break; ### not MXE
                if goodMXE: ### it's okay to write out
                  #key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(mxe[0]-1)+':'+str(mxe[1])+':'+str(f1[0]-1)+':'+str(f1[1])+':'+str(f2[0]-1)+':'+str(f2[1]);
                  key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(mxe[0]-1)+':'+str(mxe[1])+':'+str(f1[1])+':'+str(f2[0]-1);
                  if key in mxEvents: ## duplicate, 
                    dupMXE += 1;
                  else:
                    numMXEvents += 1;
                    mxEvents[key] = 1; 
                    oFile_mxe.write(str(numMXEvents)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(mxe[0]-1)+'\t'+str(mxe[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(f2[0]-1)+'\t'+str(f2[1])+'\n');
                    if (tID.find("novel")>=0): #It involves novel junction
                      neFile_mxe.write(str(numMXEvents)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(mxe[0]-1)+'\t'+str(mxe[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(f2[0]-1)+'\t'+str(f2[1])+'\n');


      #### alt-3 and alt-5 events ###
      for i in range(0,len(uf)-1): ### going through the upstream flanking exons
        e=uf[i]; 
        for j in range(i+1, len(uf)):
          u=uf[j];
          if e[0]==u[0]: ## it is alt SS event, because uf is derived from SET
            #key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(e[0]-1)+':'+str(max(e[1],u[1]))+':'+str(e[0]-1)+':'+str(min(e[1],u[1]));
            key = supInfo[1]+':'+str(min(e[1],u[1]))+':'+str(max(e[1],u[1]))+':'+str(ce[0]-1);
            if supInfo[2]=='+': ## positive strand. alt-5 event
              if key in ss5Events: ## duplicate
                dupSS5 += 1;
              else:
                ss5Events[key]=1;
                num5 += 1;
                oFile_5.write(str(num5)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(e[0]-1)+'\t'+str(max(e[1],u[1]))+'\t'+str(e[0]-1)+'\t'+str(min(e[1],u[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
                if (tID.find("novel")>=0): #It involves novel junction
                  neFile_5.write(str(num5)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(e[0]-1)+'\t'+str(max(e[1],u[1]))+'\t'+str(e[0]-1)+'\t'+str(min(e[1],u[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');

            else: ## neg strand. alt-3 event
              if key in ss3Events: ## duplicate
                dupSS3 += 1;
              else:
                ss3Events[key]=1;
                num3 += 1;
                oFile_3.write(str(num3)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(e[0]-1)+'\t'+str(max(e[1],u[1]))+'\t'+str(e[0]-1)+'\t'+str(min(e[1],u[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
                if (tID.find("novel")>=0): #It involves novel junction
                  neFile_3.write(str(num3)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(e[0]-1)+'\t'+str(max(e[1],u[1]))+'\t'+str(e[0]-1)+'\t'+str(min(e[1],u[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');


      for i in range(0,len(df)-1): ### going through the downstream flanking exons
        e=df[i];
        for j in range(i+1, len(df)):
          d=df[j];
          if e[1]==d[1]: ## it is alt SS event, because uf is derived from SET
            key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(min(e[0],d[0])-1)+':'+str(e[1])+':'+str(max(e[0],d[0])-1)+':'+str(e[1]);
            key = supInfo[1]+':'+str(ce[1])+':'+str(min(e[0],d[0])-1)+':'+str(max(e[0],d[0])-1);
            if supInfo[2]=='+': ## positive strand. alt-3 event
              if key in ss3Events: ## duplicate
                dupSS3 += 1;
              else:
                ss3Events[key]=1;
                num3 += 1;
                oFile_3.write(str(num3)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(min(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(max(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
                if (tID.find("novel")>=0): #It involves novel junction
                  neFile_3.write(str(num3)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(min(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(max(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');

            else: ## neg strand. alt-5 event
              if key in ss5Events: ## duplicate
                dupSS5 += 1;
              else:
                ss5Events[key]=1;
                num5 += 1;
                oFile_5.write(str(num5)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(min(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(max(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
                if (tID.find("novel")>=0): #It involves novel junction
                  neFile_5.write(str(num5)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(min(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(max(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');

      ### Alternative First Exon and Alternative Last Exon

###      for tID in genes[gID]: ## examine each transcript to see if the given exon is in a given transcript (index should be 1 or len-2)
###        tLen = len(genes[gID][tID]); ## length of a transcript
###        if tLen<2: ## not enough exons, skip this
###          continue; ## process next transcript
###
###        if [ce[0],ce[1]] == genes[gID][tID][1]: ## current exon is the 2nd exon in the tx
###          fEx = genes[gID][tID][0]; ## firstExon, find other tx with different fEX (non overlapping)
###          if fullyContainedInInternalExon(fEx, gID): ## fEx is fully contained.. process next transcript
###            continue;
###          for ctxID in genes[gID]: ## need to examine each tx, candidate transcript id
###            if len(genes[gID][ctxID])<2: ## not enough exons, skip this
###              continue; ## process next candidate transcript
###            if [ce[0],ce[1]] == genes[gID][ctxID][1]: ## the target exon is the 2nd exon in the ctx
###              cfEx = genes[gID][ctxID][0]; ## candidate first exon
###              if cfEx[0]>fEx[1] or fEx[0]>cfEx[1] : ### non-overlapping exon with smaller coord
###
###                #### should not fully contained in an internal exon of other transcript
###                if fullyContainedInInternalExon(cfEx, gID): ## cfEx is fully contained.. process next candidate transcript
###                  continue;
###  
###                if supInfo[2]=='+': ## positive strand. AFE, alt first exon event
###                  key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':';
###                  key += str(min(fEx[0],cfEx[0])-1)+':'+str(min(fEx[1],cfEx[1]))+':';
###                  key += str(max(fEx[0],cfEx[0])-1)+':'+str(max(fEx[1],cfEx[1]));
###                  if key in afeEvents: ## already have this one..
###                    pass; ## do nothing
###                  else: ## new AFE
###                    afeEvents[key] =1;
###                    numAFE += 1;
###                    oFile_afe.write(str(numAFE)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(min(fEx[0],cfEx[0])-1)+'\t'+str(min(fEx[1],cfEx[1]))+'\t'+str(max(fEx[0],cfEx[0])-1)+'\t'+str(max(fEx[1],cfEx[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
###                else: ## neg strand. ALE, alt last exon event
###                  key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':';
###                  key += str(min(fEx[0],cfEx[0])-1)+':'+str(min(fEx[1],cfEx[1]))+':';
###                  key += str(max(fEx[0],cfEx[0])-1)+':'+str(max(fEx[1],cfEx[1]));
###                  if key in aleEvents: ## already have this one..
###                    pass; ## do nothing
###                  else: ## new ALE
###                    aleEvents[key] =1;
###                    numALE += 1;
###                    oFile_ale.write(str(numALE)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(min(fEx[0],cfEx[0])-1)+'\t'+str(min(fEx[1],cfEx[1]))+'\t'+str(max(fEx[0],cfEx[0])-1)+'\t'+str(max(fEx[1],cfEx[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
###
###
###        if [ce[0],ce[1]] == genes[gID][tID][-2]: ## current exon is the 2nd to the last exon in the tx
###          fEx = genes[gID][tID][-1]; ## lastExon, find other tx with different fEX (non overlapping)
###          if fullyContainedInInternalExon(fEx, gID): ## fEx is fully contained.. process next transcript
###            continue;
###          for ctxID in genes[gID]: ## need to examine each tx, candidate transcript id
###            if len(genes[gID][ctxID])<2: ## not enough exons, skip this
###              continue; ## process next candidate transcript
###            if [ce[0],ce[1]] == genes[gID][ctxID][-2]: ## the target exon is the 2nd exon in the ctx
###              cfEx = genes[gID][ctxID][-1]; ## candidate last exon
###              if cfEx[0]>fEx[1] or fEx[0]>cfEx[1] : ### non-overlapping exon with smaller coord
###
###                #### should not fully contained in an internal exon of other transcript
###                if fullyContainedInInternalExon(cfEx, gID): ## cfEx is fully contained.. process next candidate transcript
###                  continue;
###
###                if supInfo[2]=='-': ## negative strand. AFE, alt first exon event
###                  key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':';
###                  key += str(min(fEx[0],cfEx[0])-1)+':'+str(min(fEx[1],cfEx[1]))+':';
###                  key += str(max(fEx[0],cfEx[0])-1)+':'+str(max(fEx[1],cfEx[1]));
###                  if key in afeEvents: ## already have this one..
###                    pass; ## do nothing
###                  else: ## new AFE
###                    afeEvents[key] =1;
###                    numAFE += 1;
###                    oFile_afe.write(str(numAFE)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(max(fEx[0],cfEx[0])-1)+'\t'+str(max(fEx[1],cfEx[1]))+'\t'+str(min(fEx[0],cfEx[0])-1)+'\t'+str(min(fEx[1],cfEx[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
###                else: ## pos strand. ALE, alt last exon event
###                  key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':';
###                  key += str(min(fEx[0],cfEx[0])-1)+':'+str(min(fEx[1],cfEx[1]))+':';
###                  key += str(max(fEx[0],cfEx[0])-1)+':'+str(max(fEx[1],cfEx[1]));
###                  if key in aleEvents: ## already have this one..
###                    pass; ## do nothing
###                  else: ## new ALE
###                    aleEvents[key] =1;
###                    numALE += 1;
###                    oFile_ale.write(str(numALE)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(max(fEx[0],cfEx[0])-1)+'\t'+str(max(fEx[1],cfEx[1]))+'\t'+str(min(fEx[0],cfEx[0])-1)+'\t'+str(min(fEx[1],cfEx[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');


      ### Retained Intron events
      for i in range(0,len(uf)): ### going through the upstream flanking exons
        f1 = uf[i]; ## first flanking exon
        for tID in genes[gID]: ## examine each transcript
          if [f1[0],ce[1]] in genes[gID][tID]: ## there is an exon starts from f1 ends at ce, it is retained intron
            #key=supInfo[1]+':'+str(f1[0]-1)+':'+str(ce[1])+':'+str(f1[0]-1)+':'+str(f1[1])+':'+str(ce[0]-1)+':'+str(ce[1]);
            key=supInfo[1]+':'+str(f1[1])+':'+str(ce[0]-1);
            if key in riEvents: ## already have this skipping events
              dupRI +=1;
              continue; ## next transcript
            else: ## new key, write it
              riEvents[key]=1;
              numRI += 1;
              oFile_ri.write(str(numRI)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(f1[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
              if (tID.find("novel")>=0): #It involves novel junction
                neFile_ri.write(str(numRI)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(f1[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');


      for i in range(0,len(df)): ### going through the downstream flanking exons
        f1 = df[i]; ## first flanking exon
        for tID in genes[gID]: ## examine each transcript
          if [ce[0],f1[1]] in genes[gID][tID]: ## there is an exon starts from ce ends at f1, it is retained intron
            key=supInfo[1]+':'+str(ce[1])+':'+str(f1[0]-1);
            if key in riEvents: ## already have this skipping events
              dupRI +=1;
              continue; ## next transcript
            else: ## new key, write it
              riEvents[key]=1;
              numRI += 1;
              oFile_ri.write(str(numRI)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(f1[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\n');
              if (tID.find("novel")>=0): #It involves novel junction
                neFile_ri.write(str(numRI)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(f1[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\n');

  ### end of else: end of merging genes
### end for gID in genes   
#
logging.debug("Number of transcripts from novel junctions: %d" % numJ);
logging.debug("Number of transcripts from novel exon, left anchored: %d" % numELA);
logging.debug("Number of transcripts from novel exon, right anchored: %d" % numERA);
logging.debug("Number of transcripts from novel exon, no anchor: %d" % numEB);
#
logging.debug("Done processing each gene from dictionary to compile AS events");
#
logging.debug("Found %d exon skipping events, from dic %d" % (numSkippingEvents, len(sEvents)));
logging.debug("duplicate skipping events: %d" % dupSE);
#
logging.debug("Found %d exon MX events, from dic %d" % (numMXEvents, len(mxEvents)));
logging.debug("duplicate MXE events: %d" % dupMXE);
#logging.debug("Filtered MXE events: %d" % filteredMXE);
#
logging.debug("Found %d alt SS events" % (len(ss3Events)+len(ss5Events)));
logging.debug("There are %d alt 5 SS events and %d alt 3 SS events." % (num5,num3));
logging.debug("duplicate alt-5 SS events: %d" % dupSS5);
logging.debug("duplicate alt-3 SS events: %d" % dupSS3);
#
#
logging.debug("Found %d AFE events, from dic %d" % (numAFE, len(afeEvents)));
#logging.debug("duplicate AFE events: %d" % dupAFE);
#
logging.debug("Found %d ALE events, from dic %d" % (numALE, len(aleEvents)));
#logging.debug("duplicate ALE events: %d" % dupALE);
#
logging.debug("Found %d RI events, from dic %d" % (numRI, len(riEvents)));
logging.debug("duplicate RI events: %d" % dupRI);
#

iFile.close();
oFile_3.close();
oFile_5.close();
oFile_ce.close();
oFile_mxe.close();
#oFile_mxe_filtered.close();
#oFile_afe.close();
#oFile_ale.close();
oFile_ri.close();
#
neFile_3.close();
neFile_5.close();
neFile_ce.close();
neFile_mxe.close();
#neFile_afe.close();
#neFile_ale.close();
neFile_ri.close();
#
#############
## calculate total running time
#############
logging.debug("Program ended");
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));

sys.exit(0);
