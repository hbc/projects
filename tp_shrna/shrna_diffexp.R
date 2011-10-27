#!/usr/bin/env RScript

# Calculate differential expression of shortRNA experimental data using DEseq.
# Usage:
#   shrna_diffexp.R <in_file> <out_base>

library(plyr)
library(reshape2)
library(DESeq)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
out.base <- args[2]

# Prepare input file in data.frame suitable for DEseq
# Does 3 transformations:
#  1. Collapse multiple probes targetting the same acession by summing
#  2. Reshape the data.frame so columns correspond to replicates and conditions
#  3. Removes any samples with missing values.
# This is specific for input data file and would need to be generalized for
# other inputs.
prepareInputs <- function(infile) {
  in.data <- read.csv(infile, header=TRUE)
  #in.data <- head(in.data, 90)
  names(in.data)[2:4] = c("counts.d.3", "counts.w.3", "accession")
  names(in.data) <- tolower(names(in.data))
  #print(head(in.data))
  collapsed.data <- ddply(in.data, .(accession, replicate), function (df)
                        data.frame(sum.d.3 = sum(df$counts.d.3),
                                   sum.w.3 = sum(df$counts.w.3),
                                   gene.symbol=df$gene.symbol[1]))
  #print(head(collapsed.data))
  collapsed.data$gene.symbol <- NULL
  melt.data <- melt(collapsed.data, id=c("accession", "replicate"),
                    measured=c("sum.d.3", "sum.w.3"))
  reshape.data <- dcast(melt.data, accession ~ variable + replicate)
  row.names(reshape.data) <- reshape.data$accession
  reshape.data$accession <- NULL
  reshape.data <- na.exclude(reshape.data)
  print(head(reshape.data))
  list(counts = reshape.data,
       conditions = c(rep("d.3", ncol(reshape.data) / 2), rep("w.3", ncol(reshape.data ) / 2)))
}

# Plot the variance diagnostic plot for a condition, checking that
# the variance function models the actual variation in the data.
plotVarianceDiagnostic <- function(cds, condition, out.base) {
  diag <- varianceFitDiagnostics(cds, condition)
  diag.file <- paste(paste(out.base, "vardiagnostic", condition, sep="-"),
                     "png", sep=".")
  png(file=diag.file)
  smoothScatter(log10(diag$baseMean), log10(diag$baseVar))
  lines(log10(fittedBaseVar) ~ log10(baseMean), diag[order(diag$baseMean),], col="red")
  dev.off()
}

# Estimate variance and assess reliability of DESeq assumptions
# This will need to be updated with next Bioconductor release
estimateVariance <- function(in.data, out.base) {
  cds <- newCountDataSet(in.data$counts, factor(in.data$conditions))
  cds <- estimateSizeFactors(cds)
  # switch to this with next release
  # cds <- estimateDispersions(cds)
  cds <- estimateVarianceFunctions(cds)
  
  scv.file <- paste(out.base, "scvplot.png", sep="-")
  png(file=scv.file)
  scvPlot(cds, ylim=c(0,2))
  dev.off()
  plotVarianceDiagnostic(cds, in.data$conditions[[1]], out.base)
  plotVarianceDiagnostic(cds, in.data$conditions[[ncol(in.data$counts)]], out.base)
  cds
}

# Calculate differential expression with DEseq, producing diagnostic plots and
# CSV output file with fold change and p-values
callDifferentialExpression <- function(cds, in.data, out.base) {
  fdr.thresh <- 0.1
  res <- nbinomTest(cds, in.data$conditions[[1]], in.data$conditions[[ncol(in.data$counts)]])
  res.sig <- res[res$padj < fdr.thresh,]
  res.sig <- res.sig[order(res.sig$padj),]
  print(head(res.sig[order(res.sig$pval),]))
  res.sig <- res.sig[,c("id", "baseMeanA", "baseMeanB", "foldChange", "padj")]
  names(res.sig) <- c("accession", in.data$conditions[[1]],
                      in.data$conditions[[ncol(in.data$counts)]],
                      "foldChange", "pval")
  out.file <- paste(out.base, "diffexp.tsv", sep="-")
  write.table(res.sig, file=out.file, row.names=FALSE, sep="\t")

  mva.file <- paste(out.base, "mvaplot.png", sep="-")
  png(file=mva.file)
  plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex=0.1,
       col = ifelse(res$padj < fdr.thresh, "red", "black"))
  dev.off()
  res
}

in.data <- prepareInputs(infile)
cds <- estimateVariance(in.data, out.base)
res <- callDifferentialExpression(cds, in.data, out.base)
