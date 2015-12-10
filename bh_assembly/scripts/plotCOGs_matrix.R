install.packages("pheatmap")
library(pheatmap)

setwd("/Volumes/ody//hanage_lab/2014_PNSP/140410_PNSP_assembly_all") # directory with the results (change to fit)

# load COG presence absence data
data <- read.table("presence_absence_matrix.out", header=T, row.names=1)

# load conversion table
strains <- read.table("strain.info")
names(strains) <- c("id", "file")
strains$seqid <- paste("seq", strains$id, sep="")

# are any COGs present in all samples? TRUE=YES
any(apply(data,2, function(x) all(x==1)))

# maximum number of shared COGs
maxsharedCOGS <- max(apply(data,2, function(x) length(which(x==1))))

# subset data to most commonly shared COGS
data.sub <- data[,apply(data,2, function(x) length(which(x==1)))>(0.75*maxsharedCOGS)] # only show the most common COGs (change to fit)

# label samples with more informative names (will be filenames of samples)
row.names(data.sub) <- sub(".fastq.gz", "", strains$file[match( row.names(data.sub),as.character(strains$seqid))])

pheatmap(data.sub, color=c("white", "grey"))
