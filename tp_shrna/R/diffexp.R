library(DESeq)

#' Plot the variance diagnostic plot for a condition, checking that
#' the variance function models the actual variation in the data.
#' @export
#' @imports DESeq
plotVarianceDiagnostic <- function(cds, condition, out.base) {
  diag <- varianceFitDiagnostics(cds, condition)
  diag.file <- paste(paste(out.base, "vardiagnostic", condition, sep="-"),
                     "png", sep=".")
  png(file=diag.file)
  smoothScatter(log10(diag$baseMean), log10(diag$baseVar))
  lines(log10(fittedBaseVar) ~ log10(baseMean), diag[order(diag$baseMean),], col="red")
  dev.off()
}

#' Estimate variance and assess reliability of DESeq assumptions
#' This will need to be updated with next Bioconductor release
#' @export
#' @imports DESeq
estimateVariance <- function(in.data, out.base) {
  cds <- newCountDataSet(in.data$counts, factor(in.data$model$conditions))
  cds <- estimateSizeFactors(cds)
  # switch to this with next release
  # cds <- estimateDispersions(cds)
  cds <- estimateVarianceFunctions(cds)
  
  scv.file <- paste(out.base, "scvplot.png", sep="-")
  png(file=scv.file)
  scvPlot(cds, ylim=c(0,2))
  dev.off()
  plotVarianceDiagnostic(cds, in.data$model$conditions[[1]], out.base)
  plotVarianceDiagnostic(cds, in.data$model$conditions[[ncol(in.data$counts)]], out.base)
  cds
}

#' Calculate differential expression with DEseq, producing diagnostic plots and
#' CSV output file with fold change and p-values
#' @export
#' @imports DESeq
callDifferentialExpression <- function(cds, in.data, out.base, config) {
  res <- nbinomTest(cds, in.data$model$conditions[[1]], in.data$model$conditions[[ncol(in.data$counts)]])
  res.sig <- res[res$padj < config$fdr_thresh,]
  res.sig <- res.sig[order(res.sig$padj),]
  print(head(res.sig[order(res.sig$pval),]))
  res.sig <- res.sig[,c("id", "baseMeanA", "baseMeanB", "foldChange", "padj")]
  names(res.sig) <- c(config$id_name, in.data$model$conditions[[1]],
                      in.data$model$conditions[[ncol(in.data$counts)]],
                      "foldChange", "pval")

  mva.file <- paste(out.base, "mvaplot.png", sep="-")
  png(file=mva.file)
  plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex=0.1,
       col = ifelse(res$padj < config$fdr_thresh, "red", "black"))
  dev.off()
  res.sig
}
