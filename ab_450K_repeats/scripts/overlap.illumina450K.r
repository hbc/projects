################################################################
## LIBRARIES ##
library(RMySQL)
################################################################
##FUNCTIONS
# Function to make it UCSC MySQL database easier to query 
query <- function(...) dbGetQuery(mychannel, ...)

## Function to call BEDtools from within R, path to BEDtools binaries must be in .Rprofile
bedTools.2in<-function(functionstring,bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring, opt.string, "-a",a.file,"-b",b.file,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}
################################################################
## ENVIRONMENT VARIABLES ##
# get base directory for analysis
basedirectory=getwd()
################################################################
## DATALOAD ##
## LOAD Tim Treche's new annotation data (v3) (from https://github.com/ttriche/hm450probes)# Tim Treche's annotation as understood by John Hutchinson - Feb 22/2012
# 1) the target sequence can be either the plus or minus strand, and has no relation to the value in the "strand" column, still useful for getting coordinates
# 3) the "start" and "end" columns give the coordinates and are for the plus strand 
# 4) the "strand" column tells you which strand the probe binds to after conversion i.e. it is the C on this strand that is actually assayed for methylation
# 5) the "cpg.GRCh37" column gives the location of the C of the CpG on the PLUS strand, even if the actual C assayed was on the minus strand (it doesn't matter though, as they are complementary)
#               *
#      10  + 5' CG----  1       1  -5' CG----  10
#          - 3' GC----             +3' GC----
#                                       *
# all this means that if the "strand" annotation is a plus, the "CpG.GRCh37" location will be at the low end of the probe, i.e. the "start" of the probe co-ordinates
# if the "strand" annotation is a minus, the "CpG.GRCh37" location will be one bp away from the "end" of the probe
load("/Users/johnhutchinson/Work/resource/Illumina450K/ttriche-hm450probes-30104c5/probes.450k.v3.rda")
ttv3.450k.annots=probes.450k
rm(probes.450k)

## extract Illumina probe locations into BED file format
probes.450K.bed=data.frame(row.names=seq(1,nrow(ttv3.450k.annots)))
probes.450K.bed$chrom=ttv3.450k.annots[,"chr"]
probes.450K.bed$chrom=paste("chr", probes.450K.bed$chrom, sep="")
probes.450K.bed$chromStart=as.numeric(ttv3.450k.annots[,"start"])
probes.450K.bed$chromEnd=as.numeric(ttv3.450k.annots[,"end"])
probes.450K.bed$name=ttv3.450k.annots[,"Probe_ID"]
probes.450K.bed$probe.type=ifelse(ttv3.450k.annots[,"design"]=="I", 1, ifelse(ttv3.450k.annots[,"design"]=="II",2, NA))
probes.450K.bed$strand=ttv3.450k.annots[,"strand"]
write.table(probes.450K.bed, file="./data/Illumina.probes.450K.ttv3.bed", quote=F, sep="\t", col.names=F, row.names=F)

############################################################################################################
## Repeat overlap
# Load repeat masker data from UCSC genome database
# using the public MySQL server for the UCSC genome browser (no password)
mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu")
# Get the UCSC repeat info, chr, start and end sites all repeats from repeatmasker track
repeat.pos.info=query("SELECT genoName, genoStart, genoEnd, repName, repClass, repFamily FROM hg19.rmsk")
# paste together repeat info so it can go through BEDtools intact
repeat.pos.info$name=paste(repeat.pos.info$repName, repeat.pos.info$repClass, repeat.pos.info$repFamily, sep="//")
repeat.pos.info.bed=repeat.pos.info[,c(1,2,3,7)]
## Use BEDtools to intersect Illumina 450K probes with repeat masker repeats and get amount of overlap
intersect.450K.repeats=bedTools.2in("intersectBED", probes.450K.bed, repeat.pos.info.bed, "-wo")
names(intersect.450K.repeats)=c("chr", "Ill.450K_start", "Ill.450K_end", "Ill.probeID", "Ill.probe.design","chr.2", "strand", "repeat_start", "repeat_end", "repeat_name//repeat_class//repeat_family", "num_bp_overlap")
# cleanup data, dropping repeat columns
intersect.450K.repeats=intersect.450K.repeats[,-grep("chr.2",names(intersect.450K.repeats))]
intersect.450K.repeats=intersect.450K.repeats[,-6]
# takeout repeat info, splitup into a format that is easier to work with and add back 
repeat.info=as.data.frame(do.call(rbind, strsplit(as.character(intersect.450K.repeats[,grep("repeat_name",names(intersect.450K.repeats))]), "\\/\\/")))
names(repeat.info)=c("repeat_name","repeat_class","repeat_family")
intersect.450K.repeats=cbind(intersect.450K.repeats, repeat.info)
intersect.450K.repeats=intersect.450K.repeats[,-grep("\\/",names(intersect.450K.repeats))]
# output to file
write.table(intersect.450K.repeats, file="./results/intersect.450K.repeats.xls", quote=F, sep="\t", col.names=T, row.names=F)

#####################################################################################################################
## SNP overlap
# Get the GATK SNP vcf files for dbSNP-excluded, 1000genomes and HapMap SNPs
download.file("ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/1.2/hg19/1000G_omni2.5.hg19.sites.vcf.gz", destfile="./data/1000G_omni2.5.hg19.sites.vcf.gz")
download.file("ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/1.2/hg19/1000G_omni2.5.hg19.vcf.gz", destfile="./data/1000G_omni2.5.hg19.vcf.gz")
download.file("ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/1.2/hg19/dbsnp_132.hg19.excluding_sites_after_129.vcf.gz", destfile="./data/dbsnp_132.hg19.excluding_sites_after_129.vcf.gz")
download.file("ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/1.2/hg19/hapmap_3.3.hg19.sites.vcf.gz", destfile="./data/hapmap_3.3.hg19.sites.vcf.gz")
system("gzip -d ./data/*.gz")


##load in 1000Genomes allele frequencies
setwd("./data")
system("sed '/##/d' 1000G_omni2.5.hg19.sites.vcf | sed 's/#//g' >1000G_omni2.5.hg19.sites.tab")
system("vcftools --vcf 1000G_omni2.5.hg19.vcf --freq --out 1000G_omni2.5.hg19")
system("paste 1000G_omni2.5.hg19.sites.tab 1000G_omni2.5.hg19.frq >1000G_omni2.5.hg19.sites.frq.tab")
system("cut -f1,2,3,4,5,11,13,14 1000G_omni2.5.hg19.sites.frq.tab >1000G_omni2.5.hg19.sites.frq.tab.1")
system("gsed 's/\t/\/\//7' 1000G_omni2.5.hg19.sites.frq.tab.1 >1000G_omni2.5.hg19.sites.frq.tab.2")
system("mv 1000G_omni2.5.hg19.sites.frq.tab.2 1000G_omni2.5.hg19.sites.frq.tab")
system("rm 1000G_omni2.5.hg19.sites.frq.tab.1")
setwd("../")
g1000.sites.freqs <- read.table("./data/1000G_omni2.5.hg19.sites.frq.tab", header=T)

#convert to BED format
g1000.sites.freqs.bed <- cbind(g1000.sites.freqs[,1:2]) 
g1000.sites.freqs.bed$END <- g1000.sites.freqs$POS+1
g1000.sites.freqs.bed$INFO <- paste(g1000.sites.freqs$ID,g1000.sites.freqs$REF,g1000.sites.freqs$ALT, g1000.sites.freqs$X.ALLELE.FREQ., sep="//")
# use bedtools function to intersect the 1000 genomes data with the Illumina data
g1000.450K.intersect <- bedTools.2in("intersectBED", probes.450K.bed, g1000.sites.freqs.bed, "-wo")

# cleanup file and extract SNPinformation
names(g1000.450K.intersect) <- c("chr","Ill.450K.low","Ill.450K.high","Ill.probeID","Ill.probe.design", "strand", "CHROM", "POS", "END", "SNPinfo", "remove")
g1000.450K.intersect <- g1000.450K.intersect[,-(grep("remove|CHROM", names(g1000.450K.intersect)))]
SNP.info=as.data.frame(do.call(rbind, strsplit(as.character(g1000.450K.intersect[,grep("SNPinfo",names(g1000.450K.intersect))]), "\\/\\/")))
names(SNP.info)=c("SNPid","REF","ALT", "REF.freq", "ALT.freq")
g1000.450K.intersect=cbind(g1000.450K.intersect, SNP.info)
g1000.450K.intersect <- g1000.450K.intersect[,-(grep("END|SNPinfo", names(g1000.450K.intersect)))]
g1000.450K.intersect$REF.freq <- do.call(rbind, strsplit(as.character(unlist(g1000.450K.intersect$REF.freq)), ":"))[,2]
g1000.450K.intersect$ALT.freq <- do.call(rbind, strsplit(as.character(unlist(g1000.450K.intersect$ALT.freq)), ":"))[,2]

# use BEDtools to intersect Illumina probe coords with SNP coords for HapMap and dbSNP vcf files (which already have allele freq info) (run system call out of R and load files back in for processing)
system("intersectBED -wo -a ./data/Illumina.probes.450K.ttv3.bed -b ./data/hapmap_3.3.hg19.sites.vcf > ./data/hapmap.tempout")
hapmap.450K.intersect=read.delim("./data/hapmap.tempout", header=F)
names(hapmap.450K.intersect) <- c("chr","Ill.450K.low","Ill.450K.high","Ill.probeID","Ill.probe.design", "strand", "CHROM", "POS","SNPid","REF", "ALT","remove1", "remove2",  "SNPinfo", "remove3")
hapmap.450K.intersect <- hapmap.450K.intersect[,-(grep("remove|CHROM", names(hapmap.450K.intersect)))]

system("intersectBED -wo -a ./data/Illumina.probes.450K.ttv3.bed -b ./data/dbsnp_132.hg19.excluding_sites_after_129.vcf > ./data/dbSNP.tempout")
dbSNP.450K.intersect=read.delim("./data/dbSNP.tempout", header=F)
names(dbSNP.450K.intersect) <- c("chr","Ill.450K.low","Ill.450K.high","Ill.probeID","Ill.probe.design", "strand", "CHROM", "POS","SNPid","REF", "ALT","remove1", "remove2",  "SNPinfo", "remove3")
dbSNP.450K.intersect <- dbSNP.450K.intersect[,-(grep("remove|CHROM", names(dbSNP.450K.intersect)))]

# process BEDtools results to determine how many bp SNP is from the CpG (adjusting for +/- strand annotation of Illumina probe - see above)
for(SNPtype in c("g1000","dbSNP","hapmap")) {
  data.temp=get(paste(SNPtype, "450K.intersect", sep=".")) 
  assign(paste(SNPtype, "distance.to.CpG",sep="."), apply(data.temp, 1, function(n) {
    strand=n[6]
    ill.start=n[2]
    ill.end=n[3]
    snp.pos=n[7]
    if (strand=="+") {
      zero=ill.start
    } else if (strand=="-") {
      zero=ill.end
    } else {
      zero=ill.end
    }
    distance=abs(as.numeric(zero)-as.numeric(snp.pos))
    return(distance)
  }))
}
# bind the distance data to the intersection data
dbSNP.450K.intersect$SNP.distance.from.CpG=dbSNP.distance.to.CpG
hapmap.450K.intersect$SNP.distance.from.CpG=hapmap.distance.to.CpG
g1000.450K.intersect$SNP.distance.from.CpG=g1000.distance.to.CpG

# parse SNPinfo column for dbSNP and hapmap dataframes to get allele frequencies (if present)
for(SNPtype in c("dbSNP","hapmap")) {
  data.temp=get(paste(SNPtype, "450K.intersect", sep="."))
  
  if (SNPtype=="dbSNP") {
    data.temp$G5=grepl("G5;", data.temp$SNPinfo)
    data.temp$G5A=grepl("G5A;", data.temp$SNPinfo)
    data.temp$ALT_freq=apply(data.temp, 1, function(y) {
      if (grepl("GMAF=", y[11])==FALSE ) {
        NA
      } else {
        info.temp.split=strsplit(y[11], ";")
        as.character(sub("GMAF=", "", info.temp.split$SNPinfo[grep("GMAF", info.temp.split$SNPinfo)]))
      }
    }) 
  } else if (SNPtype=="hapmap") {
    data.temp$ALT_freq=apply(data.temp, 1, function(x) {
      if(grepl("AF=", x[11])==FALSE ) {
        NA
      } else {
        info.temp.split=strsplit(x[11], ";")
        as.character(sub("AF=", "", info.temp.split$SNPinfo[grep("AF", info.temp.split$SNPinfo)]))
      }
    })
    data.temp$pop=apply(data.temp, 1, function(x) {
      if(grepl("set=", x[11])==FALSE ) {
        NA
      } else {
        info.temp.split=strsplit(x[11], ";")
        as.character(sub("set=", "", info.temp.split$SNPinfo[grep("set", info.temp.split$SNPinfo)]))
      }
    })
  }
  assign(paste(SNPtype, "450K.intersect", sep="."), data.temp)
}
dbSNP.450K.intersect <- dbSNP.450K.intersect[,-(grep("SNPinfo", names(dbSNP.450K.intersect)))]
hapmap.450K.intersect <- hapmap.450K.intersect[,-(grep("SNPinfo", names(hapmap.450K.intersect)))]

# reorder columns
g1000.450K.intersect=g1000.450K.intersect[,c(4,1:3,5,8,7,9:ncol(g1000.450K.intersect))]
dbSNP.450K.intersect=dbSNP.450K.intersect[,c(4,1:3,5,8,7,9:10, 12:ncol(dbSNP.450K.intersect),11)]
hapmap.450K.intersect=hapmap.450K.intersect[,c(4,1:3,5,8,7,9:10,12:ncol(hapmap.450K.intersect), 11)]
# sort by chromosome and Illumina probe start position
g1000.450K.intersect=g1000.450K.intersect[order(g1000.450K.intersect$chr, as.numeric(g1000.450K.intersect$Ill.450K.low)),]

dbSNP.450K.intersect=dbSNP.450K.intersect[order(dbSNP.450K.intersect$chr, as.numeric(dbSNP.450K.intersect$Ill.450K.low)),]
hapmap.450K.intersect=hapmap.450K.intersect[order(hapmap.450K.intersect$chr, as.numeric(hapmap.450K.intersect$Ill.450K.low)),]

### population specific MAFs
## alot done outside this script, subsetted the 1000genomes vcf file by population and then got the frqs, then meerged with GATK sites.vcf file to get SNPids and then subsetted on SNPs 
POPs <- c("AFR","AMR", "ASN", "EUR","ASW","CEU","CHB","CHS","CLM","FIN","GBR","IBS","JPT","LWK","MXL","PUR","TSI","YRI")
for (POP in POPs) {
  g1000.freqs.pop <- read.table(paste("./data/freqs", POP, "tab", sep="."), header=T)
  freq.info.pop=cbind(g1000.freqs.pop$ID, as.data.frame(do.call(rbind, strsplit(as.character(g1000.freqs.pop[,grep("X.ALLELE.FREQ.",names(g1000.freqs.pop))]), "\\/\\/"))))
  names(freq.info.pop)=c("SNPid", paste(POP, c("REF", "ALT"),  "freq", sep="."))
  freq.info.pop[,paste(POP, "REF.freq", sep=".")] <- do.call(rbind, strsplit(as.character(unlist(freq.info.pop[,paste(POP, "REF.freq",sep=".")])), ":"))[,2]
  freq.info.pop[,paste(POP, "ALT.freq", sep=".")] <- do.call(rbind, strsplit(as.character(unlist(freq.info.pop[,paste(POP, "ALT.freq",sep=".")])), ":"))[,2]
  g1000.450K.intersect <- cbind(g1000.450K.intersect, freq.info.pop[match(g1000.450K.intersect$SNPid, freq.info.pop$SNPid),2:3])
}
## calculate MAF for all pops
g1000.450K.intersect$MAF <- apply(g1000.450K.intersect[,c("REF.freq", "ALT.freq")], 1, min)
for (POP in POPs){
  g1000.450K.intersect[,paste(POP, "MAF", sep=".")] <- apply(g1000.450K.intersect[,grep(POP, names(g1000.450K.intersect))], 1, min)
}



# output to files
write.table(dbSNP.450K.intersect, file="./results/Illumina.450K.intersect.dbsnp_132.hg19.excluding_sites_after_129.JH2.xls", quote=F, sep="\t", row.names=F, col.names=T)
write.table(g1000.450K.intersect, file="./results/Illumina.450K.intersect.1000G_omni2.5.hg19.JH4.xls", quote=F, sep="\t", row.names=F, col.names=T)
write.table(hapmap.450K.intersect, file="./results/Illumina.450K.intersect.hapmap_3.3.hg19.JH2.xls", quote=F, sep="\t", row.names=F, col.names=T)

pop.data <- read.table("./data/phase1_integrated_calls.20101123.ALL.panel.panel", fill=T, header=F)
pop.data.sum <- cbind(cast(melt(pop.data), V3 ~ V2), apply(cast(melt(pop.data), V3 ~ V2), 1, sum))
pop.data.sum <- pop.data.sum[-1,-2]
names(pop.data.sum)[ncol(pop.data.sum)] <- "TOTAL"
names(pop.data.sum)[1] <- "POP"
write.table(pop.data.sum, file="./results/pop.summary.xls", col.names=T, row.names=F, sep="\t", quote=F)


