library(DESeq)
in_file = "all.combined.counts"
x = read.table(in_file, header=TRUE)
design = data.frame(row.names = colnames(x), condition = c("wt", "wt", "wt", "tsc", "tsc", "tsc"), laneType = c("single", "multi", "multi", "single", "multi", "multi"))
cds = newCountDataSet(x, design)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

simple_design = c("wt", "wt", "wt", "tsc", "tsc", "tsc")
simple_cds = newCountDataSet(x, simple_design)
simple_cds = estimateSizeFactors(simple_cds)
simple_cds = estimateDispersions(simple_cds)
simple_res = nbinomTest(simple_cds, "wt", "tsc")

plotDispEsts(cds)

# this method gives some hits for IEDs, pretty good!
fit1 = fitNbinomGLMs(cds, count ~ laneType + condition)
fit0 = fitNbinomGLMs(cds, count ~ laneType)
pvalsGLM = nbinomGLMTest(fit1, fit0)
padjGLM = p.adjust(pvalsGLM, method="BH")

rs = rowSums(counts(cds))
theta = 0.4
use = (rs > quantile(rs, probs=theta))

plot(rank(rs)/length(rs), -log10(pvalsGLM), pch=16, cex=0.45)
h1 = hist(pvalsGLM[!use], breaks=50, plot=FALSE)
h2 = hist(pvalsGLM[use], breaks=50, plot=FALSE)
colori = c("do not pass" = "khaki", "pass" = "powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori,
        space = 0, main = "", ylab="frequency")

# try the simple one
#x = read.table(in_file, header=TRUE)
#design = c("wt", "wt", "wt", "tsc", "tsc", "tsc")

cdsBlind = estimateDispersions( simple_cds, method="blind" )
vsd = varianceStabilizingTransformation( cdsBlind )
library("vsn")
par(mfrow=c(1,2))
notAllZero = (rowSums(counts(cdsBlind))>0)
meanSdPlot(log2(counts(cdsBlind)[notAllZero, ] + 1), ylim = c(0,2.5))
meanSdPlot(vsd[notAllZero, ], ylim = c(0,2.5))
mod_lfc = (rowMeans( exprs(vsd)[, pData(cdsBlind)$condition=="wt", drop=FALSE]) -
           rowMeans( exprs(vsd)[, pData(cdsBlind)$condition=="tsc", drop=FALSE]))

orderInPlot = order(pvalsGLM)
showInPlot = (pvalsGLM[orderInPlot] <= 0.08)
alpha = 0.1
plot(seq(along=which(showInPlot)), pvalsFilt[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))
abline(a=0, b=alpha/length(pvalsGLM), col="red3", lwd=2)

lfc = simple_res$log2FoldChange
logdecade = 1 + round( log10( 1+rowMeans(counts(cdsBlind, normalized=TRUE)) ) )
lfccol = colorRampPalette( c( "gray", "blue" ) )(6)[logdecade]
ymax = 4.5
pdf("log_vs_moderated.pdf")
plot( pmax(-ymax, pmin(ymax, lfc)), mod_lfc,
     xlab = "ordinary log-ratio", ylab = "moderated log-ratio",
      cex=0.45, asp=1, col = lfccol,
     pch = ifelse(lfc<(-ymax), 60, ifelse(lfc>ymax, 62, 16)))
abline( a=0, b=1, col="red3")
dev.off()

vsd2 = varianceStabilizingTransformation(cds)
library("arrayQualityMetrics")
arrayQualityMetrics(expressionset=vsd2, out_dir = "simple_report",
                    force=TRUE, do.logtransform=TRUE)
arrayQualityMetrics(expressionset=vsd2, out_dir = "simple_report",
                    intgroup=c("condition", "laneType"),
                    force=TRUE, do.logtransform=TRUE)
