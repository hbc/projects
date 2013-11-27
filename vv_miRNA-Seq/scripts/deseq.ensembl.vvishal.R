## @knitr libraries
library(DESeq)
library(plyr)
library(reshape)
library(ggplot2)
library(xtable)
library(biomaRt)
library(scales)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


## @knitr variables
if (file.exists("/n/hsphS10/hsphfs1/hsci/Vaidya/vv_miRNA-Seq")) {
  baseDir <- "/n/hsphS10/hsphfs1/hsci/Vaidya/vv_miRNA-Seq"
} else {
  baseDir <- "/Volumes/ody/consults/vv_miRNA-Seq"
}
dataDir <- file.path(baseDir, "results/htseq-count/")
resultsDir <- file.path(baseDir, "results/deseq")
metaDir <- file.path(baseDir, "meta")

count.file <- file.path(dataDir, "combined.counts")
pvalcutoff=0.2
numsig=10

gene_symbol = 'mgi_symbol'
ensembl_gene = 'mmusculus_gene_ensembl'
filter_type = 'ensembl_gene_id'


## @knitr functions
annotate_df = function(d) {
  require(biomaRt)
	ensembl = useMart('ensembl', dataset = ensembl_gene)
	a = getBM(attributes=c(filter_type, gene_symbol, "description"),
		filters=c(filter_type), values=d[, 'id'],
		mart=ensembl)
	m = merge(d, a, by.x='id', by.y=filter_type)
	return(m)
}

plotDispEsts = function(cds) {
  estimates = data.frame(means = rowMeans(counts(cds, normalized=TRUE)),
		variance = fitInfo(cds)$perGeneDispEsts)
	xg = 10^seq(-0.5, 5, length.out=300)
	yg = fitInfo(cds)$dispFun(xg)
	fitline = data.frame(xg=xg, yg=yg)
	p = ggplot(estimates, aes(means, variance)) + geom_point(size=1, alpha=0.4) +
		scale_x_log10() + scale_y_log10() +
		geom_line(data=fitline, aes(xg, yg), color="red") +
		labs(title="dispersion estimation while pooling all samples") +
		xlab("mean number of mapped reads per gene") +
		ylab("estimated dispersion")
	p
}

lm_eqn = function(df){
    m = lm(rep.2 ~ rep.1, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 2), 
              b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

qq = function(pvaldf,  cutoffpval, samples) {
  title=paste("Quantile-quantile plot of p-values", samples, sep=" - ")
  pvaldf <- pvaldf[order(pvaldf$pval, decreasing=F),]
  pvals <- as.vector(unlist(pvaldf$pval))
  padjs <- as.numeric(as.vector(unlist(pvaldf$padj)))
  colors <- as.vector(ifelse(padjs<cutoffpval, "sig", "nonsig"))
  o = -log10(pvals)
  e = -log10( 1:length(o)/length(o) )
  plot=qplot(e,o, color=colors, xlim=c(0,max(e[!is.na(e)])), ylim=c(0,max(o[!is.na(o)]))) + stat_abline(intercept=0,slope=1, col="darkgrey")
  plot=plot+labs(title=title)
  plot=plot+scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
  plot=plot+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
  plot=plot + scale_colour_manual(name="BFH adjusted pvalue", values=c("black", "red"), labels=c(paste("q>", cutoffpval, sep=""),paste("q<", cutoffpval,sep=""))) 
  plot
}

plotDE <- function(res, cutoffpval, samples ) {
  title=paste("M-A plot of", samples, sep=" - ")
  res$colors <- ifelse(res$padj<cutoffpval, "sig", "nonsig" )
  plot <- ggplot(data=res, aes(x=log(baseMean), y=log2(foldChange), colour=colors)) + 
    geom_point( size=3)  +  
    scale_colour_manual(name="BFH adjusted pvalue", values=c("#00000033","#FF0000FF"),labels=c(paste("q>", cutoffpval, sep=""),paste("q<", cutoffpval,sep=""))) +
    labs(title=title)
  plot
}


## @knitr dataload_and_reshape
counts <- read.table(file.path(dataDir, "combined.counts"), header=T, row.names=1)
counts <- counts[, order(names(counts), decreasing=T)]
# drop directories and additional trailing info from sample names, leave covariate and replicate info only 
names(counts) <- sub("results.htseq.count.", "", names(counts))
names(counts) <- sub("_.+$","", names(counts) )

# use samplenames get covars from the filenames of the counted samples
# identical values in this vector will be marked as replicates
covars <- factor(sub("rep.+$", "", names(counts)))

## load up new count dataset
cds <- newCountDataSet(counts, covars)


## @knitr rawcounts, results='asis'
annotated_counts = head(as.data.frame(counts(cds)), numsig)
annotated_counts$id = rownames(annotated_counts)
annotated_counts = annotate_df(annotated_counts)
print(xtable(annotated_counts), type='html', include.rownames=F)


## @knitr sizefactors, tidy=TRUE
cds = estimateSizeFactors(cds)
sizeFactors(cds)


## @knitr normalized_counts, results='asis'
annotated_normalized = head(as.data.frame(counts(cds, normalized=TRUE)), numsig)
annotated_normalized$id = rownames(annotated_normalized)
annotated_normalized = annotate_df(annotated_normalized)
print(xtable(annotated_normalized), type='html', include.rownames=F)


## @knitr ratio_hist, fig.cap=""
raw.counts = counts(cds, normalized=FALSE)
cols <- sample(ncol(raw.counts),2, replace=F)
rawdata = data.frame(ratio=raw.counts[,cols[1]] / raw.counts[,cols[2]])
rawdata$set <- "raw"
norm.counts = counts(cds, normalized=TRUE)
normdata = data.frame(ratio=norm.counts[,1] / norm.counts[,2])
normdata$set <- "normalized"
raw.norm.data <- rbind(rawdata, normdata)

n = ggplot(raw.norm.data, aes(x=ratio, fill=set)) + geom_density(alpha=0.25) +
	scale_x_log10(breaks=c(0.01, 0.1, 1, 10, 100), labels=math_format(format=log10)) +
	labs(title="Normalized counts")
print(n)


## @knitr estimate_sizefactors, results='hide', fig.cap="Empirical and fitted dispersion values plotted against mean expression strength"
# sharingMode = maximum, most conservative approach to sharing information across genes to reduce variability of the dispersion estimates
cds <- estimateDispersions(cds, method="pooled", sharingMode="maximum", fitType="parametric")
plotDispEsts(cds)


## @knitr pca
pcaPlot <- function(readdata, title)  {
  fit <- prcomp(t(readdata))
  colors <- cbPalette[factor(pData(cds)$condition)]
  legend_values=unique(cbind(colors, as.character(pData(cds)$condition)))
  ##all samples
  plot(fit$x, bg=colors, col="black", cex=2,pch=21, main=title, oma=c(8,5,5,14))
  legend("topright", cex=0.7, col="black", pt.bg=legend_values[,1], pt.cex=1.25, legend=legend_values[,2],  pch=21, bty="n", x.intersp=1)
  }
pcaPlot(raw.counts, "Raw counts")
pcaPlot(norm.counts, "Normalized counts")


## @knitr filter
  ## get sum of counts for all samples for each gene
  rowcounts <- rowSums(norm.counts)
  ## filter the data based on the minimal row sum 
  use <- (rowcounts > quantile(rowcounts, 0.4))
  cds.filt <- cds[use,]
  filt.norm.counts <- norm.counts[use,]


## @knitr pairwise_comparisons, cache=TRUE
## first construct the actual combinations
all.pair.combos <- combn(as.vector(unique(pData(cds)$condition)),2)
vs.norm.combos <- all.pair.combos[,which(all.pair.combos[1,]=="Normal")]
setnames <- apply(vs.norm.combos, 2, function(n) paste(n[1], n[2], sep="-vs-"))
sig.results <- alply(vs.norm.combos, 2, function(combo) {
  setname <- paste(combo[1], combo[2], sep="-vs-")
  print(setname)
  ## perform significance testing
  res.filt <- nbinomTest(cds.filt, combo[1], combo[2])
  ## get normalized counts for significant hits, relabel samples with condition rather than sampleID
  filt.norm.counts.d <- as.data.frame(filt.norm.counts)
  results.1 <- filt.norm.counts.d[which(res.filt$padj<pvalcutoff),]
  ## get means and pvalues for significant hits and put together with counts
  results.2 <- res.filt[which(res.filt$padj<pvalcutoff),]
  results <- cbind(results.1, results.2)
  results <- annotate_df(results)
  results <- results[order(results$padj),]
  ## output some plots
  qqplots <- qq(res.filt[,c("pval", "padj")], pvalcutoff, setname)
  DEplots <- plotDE(res.filt, pvalcutoff, setname)
  return(list(res.filt=res.filt, results=results,  qqplots=qqplots, DEplots=DEplots))
})  


## @knitr out1, fig.width=11, fig.height=6
sig.results[[5]]$qqplots
sig.results[[5]]$DEplots


## @knitr tables1, results='asis'
if (nrow(sig.results[[5]]$results)>(numsig-1)) {
  out <- xtable(sig.results[[5]]$results[1:numsig,])
  } else {
    out  <- xtable(sig.results[[5]]$results)
    }
print(out, type='html',include.rownames=FALSE)
write.table(sig.results[[5]]$results, file=file.path(resultsDir, paste("DE.genes.q", pvalcutoff, setnames[5], "xls", sep=".")), quote=F, sep="\t", row.names=F, col.names=T)


## @knitr out2, fig.width=11, fig.height=6
sig.results[[3]]$qqplots
sig.results[[3]]$DEplots


## @knitr tables2, results='asis'
if (nrow(sig.results[[3]]$results)>(numsig-1)) {
  out <- xtable(sig.results[[3]]$results[1:numsig,])
  } else {
    out  <- xtable(sig.results[[3]]$results)
    }
print(out, type='html',include.rownames=FALSE)
write.table(sig.results[[3]]$results, file=file.path(resultsDir, paste("DE.genes.q", pvalcutoff, setnames[3], "xls", sep=".")), quote=F, sep="\t", row.names=F, col.names=T)


## @knitr out3, fig.width=11, fig.height=6
sig.results[[2]]$qqplots
sig.results[[2]]$DEplots


## @knitr tables3, results='asis'
if (nrow(sig.results[[2]]$results)>(numsig-1)) {
  out <- xtable(sig.results[[2]]$results[1:numsig,])
  } else {
    out  <- xtable(sig.results[[2]]$results)
    }
print(out, type='html',include.rownames=FALSE)
write.table(sig.results[[2]]$results, file=file.path(resultsDir, paste("DE.genes.q", pvalcutoff, setnames[2], "xls", sep=".")), quote=F, sep="\t", row.names=F, col.names=T)


## @knitr out4, fig.width=11, fig.height=6
sig.results[[1]]$qqplots
sig.results[[1]]$DEplots


## @knitr tables4, results='asis'
if (nrow(sig.results[[1]]$results)>(numsig-1)) {
  out <- xtable(sig.results[[1]]$results[1:numsig,])
  } else {
    out  <- xtable(sig.results[[1]]$results)
    }
print(out, type='html',include.rownames=FALSE)
write.table(sig.results[[1]]$results, file=file.path(resultsDir, paste("DE.genes.q", pvalcutoff, setnames[1], "xls", sep=".")), quote=F, sep="\t", row.names=F, col.names=T)


## @knitr out5, fig.width=11, fig.height=6
sig.results[[6]]$qqplots
sig.results[[6]]$DEplots


## @knitr tables5, results='asis'
if (nrow(sig.results[[6]]$results)>(numsig-1)) {
  out <- xtable(sig.results[[6]]$results[1:numsig,])
  } else {
    out  <- xtable(sig.results[[6]]$results)
    }
print(out, type='html',include.rownames=FALSE)
write.table(sig.results[[6]]$results, file=file.path(resultsDir, paste("DE.genes.q", pvalcutoff, setnames[6], "xls", sep=".")), quote=F, sep="\t", row.names=F, col.names=T)


## @knitr out6, fig.width=11, fig.height=6
sig.results[[4]]$qqplots
sig.results[[4]]$DEplots


## @knitr tables6, results='asis'
if (nrow(sig.results[[4]]$results)>(numsig-1)) {
  out <- xtable(sig.results[[4]]$results[1:numsig,])
  } else {
    out  <- xtable(sig.results[[4]]$results)
    }
print(out, type='html',include.rownames=FALSE)
write.table(sig.results[[4]]$results, file=file.path(resultsDir, paste("DE.genes.q", pvalcutoff, setnames[4], "xls", sep=".")), quote=F, sep="\t", row.names=F, col.names=T)


## @knitr save_image
save.image(file.path(resultsDir, "RDATA" ))


