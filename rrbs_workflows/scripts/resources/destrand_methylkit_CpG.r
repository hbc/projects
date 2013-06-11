# destrand output from methylKit
# drafted by Allan Just 11/8/2012

# merges in strand info from a Bismark Genome-Wide CpG report 
#  (hg19 -see bismark_methylation_extractor script)
# then uses code from Methylkit to destrand
# operates in a for loop without error catching
# stores a few summary fields in a tab delimited output "destrand_summary_YYYY-MM-DD.txt"

# here 'datadir' contains "methylation_quantitation_results/" folder
datadir <- "/n/hsphS10/hsphfs1/chb/projects/rrbs_workflows/sample_data/AB_data/no_coverage_filter/"
# 'resourcedir' contains the R object with genome-wide CpGs ~56 million cytosines
resourcedir <- "/n/hsphS10/hsphfs1/chb/projects/rrbs_workflows/scripts/resources/"

library(methylKit)
library(data.table)
library(plyr)

# function taken from methylKit
# unifies forward and reverse strand CpGs on the forward strand if both are on the same CpG
# if that's the case their values are generally correlated
CpG.dinuc.unify<-function(cpg)
{
  
  cpgR=cpg[cpg$strand=="-",]
  cpgF=cpg[cpg$strand=="+",]
  cpgR$start=cpgR$start-1
  cpgR$end=cpgR$end-1
  cpgR$strand="+"
  
  cpgR$id=paste(cpgR$chr,cpgR$start,sep=".")
  
  cpgFR=merge(cpgF,cpgR,by="id")
  #hemi =cpgFR[abs(cpgFR$freqC.x-cpgFR$freqC.y)>=50,]
  #cpgFR=cpgFR[abs(cpgFR$freqC.x-cpgFR$freqC.y)<50,]
  res=data.frame(
    id =as.character(cpgFR$id),
    chr =as.character(cpgFR$chr.x),
    start =(cpgFR$start.x),
    end =(cpgFR$start.x),
    strand =rep("+",nrow(cpgFR)),
    coverage=cpgFR$coverage.x + cpgFR$coverage.y,
    numCs =cpgFR$numCs.x + cpgFR$numCs.y ,
    numTs =cpgFR$numTs.x + cpgFR$numTs.y ,stringsAsFactors =F
  )
  res=rbind(res, cpgF[ !cpgF$id %in% res$id,],cpgR[ !cpgR$id %in% res$id,] )
  res=res[order(res$id),]
  return(res)
}

# turned off - read in the genome_wide_CpG_report from text
if(0){
  genomecpgfile <- "/n/hsphS10/hsphfs1/chb/projects/rrbs_workflows/sample_data/AB_data/B1627_GGCTAC_L002_R1.trimmed.fq_bismark.coordsorted.genome_wide_CpG_report.txt"
  system.time(allcpg <- read.table(genomecpgfile, header = F, nrows = 56434896, 
    colClasses = c("factor", "integer", "factor", 
                   rep("NULL", 4))) )
  object.size(allcpg)
  names(allcpg) <- c("chr", "site", "strand")
  save(allcpg, file = "Genome-wide_CpG.Rdata")
}

# stored as an R object for faster input
system.time(load(paste(resourcedir, "Genome-wide_CpG.Rdata", sep = "")))

# convert to data.table for faster operations
allcpg <- data.table(allcpg)
# change site to start to match methylRaw files
setnames(allcpg, "site", "start")
setkey(allcpg,chr,start)

# loop next section across all
# methylKit output
methylKitfilenames <- list.files(paste(datadir, "methylation_quantitation_results/", sep = ""), pattern = "R1_CpG\\.txt")

# store some stats
destrandlist <- list()

for(i in 1:length(methylKitfilenames)){
# import a methylkit object 
sampfile <- methylKitfilenames[i]
sampname <- sub("\\.txt", "", sampfile)
methylRawsamp <- read(paste(datadir, "methylation_quantitation_results/", sampfile, sep = ""), 
  sampname, assembly = "hg19", pipeline = "bismark")
#methylRawnames <- c("id", "chr", "start", "end", "strand", "coverage", "numCs", "numTs")
## convert to data.table
methylRawsamp.dt <- data.table(getData(methylRawsamp), key = c("chr","start"))
## remove useless strand
methylRawsamp.dt[, strand := NULL]
## merge in strand from genome-wide file (inner join)
methylRawsamp.dt <- allcpg[methylRawsamp.dt]
#methylRawsampstrand <- data.frame(methylRawsampstrand.dt)
##  push back into methylRaw object
#testout <- new("methylRaw", 
#  methylRawsampstrand[, methylRawnames],
#  sample.id = "testout", assembly = "hg19", context = "CpG", resolution = "base")
##  destrand function taken from methylkit 
sampdtdestrand <- CpG.dinuc.unify(methylRawsamp.dt)
write.table(sampdtdestrand, file = paste(datadir, "methylation_quantitation_results/",sampname, "_DS.txt", sep = ""),
  row.names = F)
print(paste("destranding", sampname))
# a few quick statistics while reading in every sample
destrandlist[[sampname]] <- 
  data.frame(cpgonbothstrands = nrow(methylRawsamp.dt) - nrow(sampdtdestrand),
  covgt1x_bs = nrow(methylRawsamp.dt),
  covgt1x = nrow(sampdtdestrand),
             covge10x_bs = sum(methylRawsamp.dt$coverage >= 10),
             covge10x = sum(sampdtdestrand$coverage >= 10),
             covge50x_bs = sum(methylRawsamp.dt$coverage >= 50),
             covge50x = sum(sampdtdestrand$coverage >= 50))

rm(methylRawsamp, methylRawsamp.dt, sampdtdestrand)
gc()
} #end of loop


# collapse summary to data.frame
destrandinfo <- ldply(destrandlist)
# save summary as csv file
destrandfilename <- paste("destrand_summary_", Sys.Date(), ".csv", sep = "")
write.csv(destrandinfo, file = paste(datadir, destrandfilename,sep = ""), row.names = F)
print(paste(destrandfilename, "saved to ", datadir))

## explore sample to confirm destranding works
#methylRawsamp.dt[ , diff := c(NA, diff(start))]
#table(methylRawsamp.dt$diff == 1, methylRawsamp.dt$strand)
#nrow(methylRawsamp.dt) - nrow(sampdtdestrand)
## no sites with negative strand or with diff of 1 left
#sampdtdestrand$diff <- c(NA, diff(sampdtdestrand$start))
#table(sampdtdestrand$diff == 1, sampdtdestrand$strand)    


#End of file