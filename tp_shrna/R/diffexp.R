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
  cds <- newCountDataSet(in.data$counts, factor(in.data$model$condition))
  cds <- estimateSizeFactors(cds)
  # Check for replicated conditions
  if (length(unique(in.data$model$condition)) == length(in.data$model$condition)) {
    # No replicates, treat all samples together for variance estimation
    cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only",
                               fitType="local")
  } else {
    # Replicates, standard variance estimation
    cds <- estimateDispersions(cds)
  }
  if (!is.null(config$multifactor) && config$multifactor) {
    in.data$model$condition <- factor(in.data$model$condition)
    cds_multi <- newCountDataSet(in.data$counts, in.data$model)
    cds_multi <- estimateSizeFactors(cds_multi)
    cds_multi <- estimateDispersions(cds_multi)
  } else {
    cds_multi <- NULL
  }

  if (!is.null(config$out_base)) {
    plotDispEsts(cds, config$out_base)
    if (!is.null(cds_multi)) {
      plotDispEsts(cds_multi, paste(config$out_base, "multi", sep="-"))
    }
  }
  list(cds=cds, cds_multi=cds_multi)
}

#' Get pvalue to select by and threshold based on configuration
#' Can use either adjusted p-value (FDR) or standard p-value
#' Standard p-value is more useful in non-replicate experiments
getPvalThresh <- function(config) {
  if (!is.null(config$pval_thresh)) {
    list(attr="pval", thresh=config$pval_thresh, head=config$top_targets)
  } else {
    list(attr="padj", thresh=config$fdr_thresh, head=config$top_targets)
  }
}

#' Prepare DESeq result by filtering on p-value and assigning columns
prepareDEResult <- function(res, in.data, config) {
  pselect <- getPvalThresh(config)
  print(head(res[order(res[pselect$attr]), ]))
  if (!is.null(pselect$head)) {
    res.sig <- head(res[order(res[pselect$attr]), ], n=pselect$head)
  } else {
    res.sig <- res[res[pselect$attr] < pselect$thresh, ]
  }
  res.sig <- res.sig[order(res.sig[pselect$attr]), ]
  res.sig <- res.sig[,c("id", "baseMeanA", "baseMeanB", "foldChange", pselect$attr)]
  names(res.sig) <- c(config$id_name, in.data$model$condition[[1]],
                      in.data$model$condition[[ncol(in.data$counts)]],
                      "foldChange", "pval")
  res.sig
}

#' MvA plot of differential expression result
mvaPlot <- function(res, out_base, config) {
  pselect <- getPvalThresh(config)
  mva.file <- paste(out_base, "mvaplot.png", sep="-")
  png(file=mva.file)
  plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex=0.1,
       col = ifelse(res[pselect$attr] < pselect$thresh, "red", "black"))
  dev.off()
}

#' Calculate differential expression with DEseq, producing diagnostic plots and
#' CSV output file with fold change and p-values
#' @export
#' @imports DESeq
callDifferentialExpression <- function(cds, cds_multi, in.data, config) {
  res <- nbinomTest(cds, in.data$model$condition[[1]],
                    in.data$model$condition[[ncol(in.data$counts)]])
  res_sig <- prepareDEResult(res, in.data, config)
  if (!is.null(cds_multi)) {
    fit1 <- fitNbinomGLMs(cds_multi, count ~ type + condition)
    fit0 <- fitNbinomGLMs(cds_multi, count ~ type)
    pvalsGLM <- nbinomGLMTest(fit1, fit0)
    padjGLM <- p.adjust(pvalsGLM, method="BH")
    res_multi <- res
    res_multi$pval <- pvalsGLM
    res_multi$padj <- padjGLM
    res_multi_sig <- prepareDEResult(res_multi, in.data, config)
  } else {
    res_multi <- NULL
    res_multi_sig <- NULL
  }

  if (!is.null(config$out_base)) {
    mvaPlot(res, config$out_base, config)
    if (!is.null(res_multi)) {
      mvaPlot(res_multi, paste(config$out_base, "multi", sep="-"), config)
    }
  }
  list(res=res_sig, res_multi=res_multi_sig)
}
