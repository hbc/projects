#!/n/sw/R-2.15.0/bin/Rscript

library(methylKit)

args=commandArgs(TRUE)
filename=args[1]
sampleID=sub("sam", "",filename)

dataDir=args[2]

read.bismark(location=filename, sample.id=sampleID, assembly="mm9", save.folder=dataDir, save.context="CpG", read.context="CpG", nolap=FALSE, mincov=10, minqual=20, phred64=FALSE)

q()