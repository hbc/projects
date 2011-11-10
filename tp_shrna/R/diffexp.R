library(DESeq)

#' Plot the variance diagnostic plot for a condition, checking that
#' the variance function models the actual variation in the data.
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

#' Plot the estimated dispersions versus raw data
#' @imports DESeq
plotDispEsts <- function(cds, out_base) {
  out_file = paste(out_base, "displot.png", sep="-")
  png(file=out_file)
  plot(rowMeans(counts(cds, normalized=TRUE)),
       fitInfo(cds)$perGeneDispEsts,
       pch = ".", log = "xy")
  xg <- 10^seq(-.5, 5, length.out=300)
  lines(xg, fitInfo(cds)$dispFun(xg), col="red")
  dev.off()
}


#' Estimate variance and assess reliability of DESeq assumptions
#' This will need to be updated with next Bioconductor release
#' @export
#' @imports DESeq
estimateVariance <- function(in.data, config) {
  if (config$multifactor) {
    cds <- newCountDataSet(in.data$counts, in.data$model)
  } else {
    cds <- newCountDataSet(in.data$counts, factor(in.data$model$conditions))
  }
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)

  if (!is.null(config$out_base)) {
    plotDispEsts(cds, config$out_base)
  }
  cds
}

#' Calculate differential expression with DEseq, producing diagnostic plots and
#' CSV output file with fold change and p-values
#' @export
#' @imports DESeq
callDifferentialExpression <- function(cds, in.data, config) {
  if (config$multifactor) {
    fit1 <- fitNbinomGLMs(cds, count ~ type + condition)
    fit0 <- fitNbinomGLMs(cds, count ~ type)
    pvalsGLM <- nbinomGLMTest(fit1, fit0)
    res <- p.adjust(pvalsGLM, metho="BH")
  } else {
    res <- nbinomTest(cds, in.data$model$conditions[[1]],
                      in.data$model$conditions[[ncol(in.data$counts)]])
  }
  res.sig <- res[res$padj < config$fdr_thresh,]
  res.sig <- res.sig[order(res.sig$padj),]
  print(head(res.sig[order(res.sig$pval),]))
  res.sig <- res.sig[,c("id", "baseMeanA", "baseMeanB", "foldChange", "padj")]
  names(res.sig) <- c(config$id_name, in.data$model$conditions[[1]],
                      in.data$model$conditions[[ncol(in.data$counts)]],
                      "foldChange", "pval")

  if (!is.null(config$out_base)) {
    mva.file <- paste(config$out_base, "mvaplot.png", sep="-")
    png(file=mva.file)
    plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex=0.1,
         col = ifelse(res$padj < config$fdr_thresh, "red", "black"))
    dev.off()
  }
  res.sig
}
