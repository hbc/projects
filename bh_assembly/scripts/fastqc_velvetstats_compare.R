library(plyr)
library(reshape2)
library(ggplot2)

setwd("/Volumes/ody/bh_assembly/Gaza/")
baseDir <- getwd()


fastqs <- list.files(pattern="fastq$")
samples <- unique(sub("_L001_R._001.fastq", "", fastqs))
fastqcdirs <- file.path(getwd(), sub(".fastq", "_fastqc",fastqs ))
fastqcdirs <- t(matrix(fastqcdirs, nrow=2))
velvetdirs <- file.path(getwd(), paste(samples, "velvet", sep="_"))
velvetlogs <- file.path(baseDir, paste(samples, "velvet_assembly.log", sep="_") )

data_files <- do.call(cbind, list(samples, fastqcdirs,velvetdirs, velvetlogs))

results <- adply(data_files, 1, function(x){
  sample <- x[1]
  print(sample)
  fastqcdir1 <- x[2]
  fastqcdir2 <- x[3]
  velvetdir <- x[4]
  velvetlogfile <- x[5]
  
  fastqdatafile1 <- file.path(fastqcdir1, "fastq_data.txt")
  fastqdatafile2 <- file.path(fastqcdir2, "fastq_data.txt")
  velvetstatsfile <- file.path(velvetdir, "stats.txt")
    
  setwd(fastqcdir1)
  tempfile <- tempfile()  
  system.call <- paste('sed -n "/Base\tMean/,/END/p" fastqc_data.txt >', tempfile, sep=" ")
  system(system.call)
  fastqcdata1 <- read.delim(tempfile)
  fastqcdata1 <- fastqcdata1[!grepl("END", fastqcdata1$X.Base),]
  fastqcdatamean1 <- mean(fastqcdata1$Mean)
  
  setwd(fastqcdir2)
  tempfile <- tempfile()  
  system.call <- paste('sed -n "/Base\tMean/,/END/p" fastqc_data.txt >', tempfile, sep=" ")
  system(system.call)
  fastqcdata2 <- read.delim(tempfile)
  fastqcdata2 <- fastqcdata2[!grepl("END", fastqcdata2$X.Base),]
  fastqcdatamean2 <- mean(fastqcdata2$Mean)
  
  tempfile <- tempfile()  
  system.call <- paste('grep "Total Sequences" fastqc_data.txt >', tempfile, sep=" ")
  system(system.call)
  numreads <- readLines(tempfile)
  numreads <- as.numeric(strsplit(numreads,"\t")[[1]][2])
    
  setwd(baseDir)
  velvetlogdata <- readLines(velvetlogfile)
  kmer <- as.numeric(strsplit(velvetlogdata[length(velvetlogdata)-1], " |=")[[1]][2])
  
  setwd(velvetdir)
  myStatsTable<-read.table(velvetstatsfile,header=TRUE)
  contigs<-rev(sort(myStatsTable$lgth+kmer-1))
  n50<-contigs[cumsum(contigs) >= sum(contigs)/2][1]
  
  return(c(sample, numreads, fastqcdatamean1,fastqcdatamean2, n50))
})

results$X1 <- NULL

names(results) <- c("sample", "numreads", "mean_quality1","mean_quality2", "n50")
results$mean_quality1 <- as.numeric(results$mean_quality1)
results$mean_quality2 <- as.numeric(results$mean_quality2)
results$n50 <- as.numeric(results$n50)
results$numreads <- as.numeric(results$numreads)

ggplot(results, aes(x=mean_quality1, y=n50, size=numreads, label=sample))+
  geom_point()+scale_size(range=c(2,100))


cla