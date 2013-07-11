library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(oligo)
library(genefilter)
library(plyr) 
library(limma)
library(sva)
library(arrayQualityMetrics)
library(pheatmap)
library(Biobase)

# CEL-file handling
baseDir <- "/Users/johnhutchinson/Work/Consults/Glimcher_bone_development_13"
dataDir <- file.path(baseDir, 'data')
resultsDir <- file.path(baseDir, 'results')

# Parse array files and normalize
setwd(baseDir)
celFiles <- list.celfiles(dataDir, pattern='LG201204', full.names=TRUE)
affy <- read.celfiles(celFiles, verbose=T)

# Set up covariates based on sample information from Anju
pDataFile <- file.path(dataDir, 'sampleInfo_secondGroup.txt')
pDataObj <- read.table(pDataFile, row.names=1, header=T, sep='\t')
all(rownames(pDataObj) == colnames(exprs(affy)))
pDataObj$Batch <- factor(c(rep(1, 4), rep(2, 4)))
pData(affy) <- pDataObj
pData(affy)

# Transcript (gene) level normalization using RMA
allArrays <- rma(affy, target='core')
allArrays

#
# Correct for batch effects
#
# Create model with Condition as factor variable
mod <- model.matrix(~as.factor(Condition), data=pDataObj)
batch <- pData(allArrays)$Batch

# Modify expression matrix
edata <- exprs(allArrays)
combat_edata <- ComBat(dat=edata,
                       batch=batch,
                       mod=mod,
                       numCovs=NULL,
                       par.prior=TRUE, 
                       prior.plots=TRUE)
allArraysBatch <- allArrays
exprs(allArraysBatch) <- combat_edata

# Retrieving NetAffx Biological Annotation
featureData(allArraysBatch) <- getNetAffx(allArraysBatch, 'transcript')
varLabels(featureData(allArraysBatch))

# Extract the 'gene assignment' annotation
annot <- pData(featureData(allArraysBatch)[, c('geneassignment')])
head(annot[!is.na(annot), ], 2)

# Generate a list of gene symbols from the gene assignment
desc <- annot[, 1]
symbols <- unlist(lapply(desc, function(x) strsplit(x, ' // ')[[1]][2]))
length(featureData(allArraysBatch)$probesetid) == length(symbols)
head(symbols[!is.na(symbols)])

#
# Grabbing gene list and order from Anju's data
#
setwd(baseDir)
sheet2 <- read.table('data/20120531 Gene Lists/B_combined_sheet2.csv',
                     sep=',', header=T)
array2 <- sheet2[, c(1, 2, 3, 6, 9)]
array2 <- array2[!duplicated(array2), ]

# All we need is the ProbeID
sheet2probes <- array2$Row.names

#
# Subset expression set by probe ID
#
arraySubset <- exprs(allArraysBatch[sheet2probes, ])
colnames(arraySubset) <- pData(allArraysBatch)$Sample
rownames(arraySubset) <- paste(array2$Row.names,
                               array2$ID_21,
                               sep='.')

# Visualization
colors <- brewer.pal(9, 'RdYlGn')
pal <- rev(colorRampPalette(colors)(50))
#range <- seq(-5, 5, 0.2)

# Define colors for the annotation
annotation <- data.frame(Var1=pData(allArraysBatch)$Condition)
rownames(annotation) <- pData(allArraysBatch)$Sample
Var1 <- brewer.pal(4, 'Spectral')
names(Var1) <- levels(pData(allArraysBatch)$Condition)
ann_colors <- list(Var1=Var1)


pheatmap(arraySubset,
         color=pal,
#         breaks=range,
         cellwidth=25,
         cellheight=12,
         scale='none',
         cluster_rows=T,
         cluster_cols=T,
         legend=T,
         fontsize_col=10,
         annotation=annotation,
         annotation_legend=T,
         annotation_colors=ann_colors,
         filename='20120703_Sheet2_sorted.pdf')

