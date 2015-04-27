if(file.exists("/n/hsphS10/hsphfs1/chb/projects/mg_spikeins/")){
  baseDir="/n/hsphS10/hsphfs1/chb/projects/mg_spikeins/"
    } else if (file.exists("/Users/johnhutchinson/projects/mg_spikeins/")){
    baseDir="/Users/johnhutchinson/projects/mg_spikeins/"
    }

dataDir <- file.path(baseDir, "data")
metaDir <- file.path(baseDir, "meta")
resultsDir <- file.path(baseDir, "results")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colorblind friendly palette
covarsfilename="covars.desc" # tab delimited file describing samples
lowintensity.percentile=0.1
mad.quantile.cutoff=0.1
pvalue.cutoff=0.25



## ----libraries_variables, echo=TRUE--------------------------------------
library(knitr)
library(affy) # for loess normalization
library(oligo) # array utilities
library(pd.mogene.2.0.st) # array layout annotation
library(mogene20sttranscriptcluster.db) # array probe to gene annotations
library(arrayQualityMetrics) # array quality control reports
library(limma) # array statistical analyses
library(pheatmap) # pretty heatmaps
library(plyr) # data format utility
library(reshape2) # data format utility
library(devtools) # install libraries from github
install_git("git://github.com/hbc/CHBUtils.git") # misc personal utilities
library(CHBUtils)
install_git("git://github.com/stephaniehicks/quantro")
library(quantro)
library(ggplot2) # pretty graphs
library(ggdendro) # for pretty dendrograms
library(RColorBrewer)


## ----functions-----------------------------------------------------------
PCAplot.eset <- function(eset=NULL, categories1=NULL, categories2=NULL,title=NULL, colorpalette=NULL, alpha=1, numcomponents=4){
  alpha <- sprintf("%x", ceiling(alpha*255))
  colorpalette <- paste(colorpalette, alpha, sep="")
  eset.core <- exprs(eset) 
  myPca.core <- prcomp(t(eset.core))
  tmpPCAData.core <- as.data.frame(myPca.core$x[,1:numcomponents])
  pd <- pData(eset)
  colors <- colorpalette[factor(as.character(unlist(pd[,categories1])))]
  numcat2s <- length(unique(as.character(unlist(pd[,categories2])))) 
  shapes <- c(21:25)[1:numcat2s][factor(as.character(unlist(pd[,categories2])))]
  pairs(tmpPCAData.core, bg=colors, col="#606060", cex=2, pch=shapes, main=title, oma=c(8,5,5,14))
  legend("right", cex=0.7, col="#606060", pt.bg=colors, pt.cex=1.5, legend=paste(pd[,categories1], pd[,categories2], sep="-"),  pch=shapes, bty="n", x.intersp=1)
}

PCAplot.sd.eset <- function(eset=NULL,  title=NULL){
  eset.core <- exprs(eset)
  myPca.core <- prcomp(t(eset.core))
  # SD of components
  sdevdf <- data.frame(cbind(as.numeric(myPca.core$sdev),c(1:length(myPca.core$sdev))))
  sdevdf$prop <-  sdevdf$X1/sum(sdevdf$X1)
  sdevdf$cum <- cumsum(sdevdf$prop)
  ggplot(sdevdf, aes(x=X2, y=prop)) + 
    geom_point(size=4, color="red") + 
    scale_x_continuous('Component') + 
    scale_y_continuous('Standard Deviation') +
    ggtitle(title) +
    geom_line(data=sdevdf, aes(x=X2, y=cum))
}


## ----dataload, results='hide'--------------------------------------------
covars <- read.table(file.path(metaDir, covarsfilename),header=TRUE, sep="\t", row.names=1) # simple tab delimited file with CEL file in first column (no heading for this column) and sample metadata (i.e. sampleID, treatment group, batch etc.) in subsequent columns

celFiles <- file.path(dataDir, row.names(covars))
affyRaw <- read.celfiles(celFiles)
pData(affyRaw) <- covars 
sampleNames(affyRaw) <- pData(affyRaw)$sampleID
validObject(affyRaw)
rm(covars)


## ----covars, results='asis'----------------------------------------------
# Sample information table
kable(pData(affyRaw))


## ----rawQC, eval=FALSE---------------------------------------------------
## arrayQualityMetrics(expressionset=affyRaw, outdir=file.path(resultsDir, 'report_raw'), force=TRUE, do.logtransform=TRUE, intgroup=c("genotype", "treatment"))


## ----normalize, results='hide'-------------------------------------------
affyNorm.core <- rma(affyRaw, target="core", background=TRUE, normalize=TRUE)


## ----normQC, eval=FALSE--------------------------------------------------
## arrayQualityMetrics(expressionset=affyNorm.core, outdir=file.path(resultsDir, paste("report_rma.core", sep=".")), force=TRUE, do.logtransform=FALSE)


## ----removeSHsamples-----------------------------------------------------
affyNorm.core <- affyNorm.core[,which(pData(affyNorm.core)$treatment!="SH")]


## ----cluster1, out.width='100%'------------------------------------------
plot_dendro(affyNorm.core, title="", labels.colname="treatment", colors.colname="genotype")


## ----PCAsd1, out.width='75%'---------------------------------------------
PCAplot.sd.eset(affyNorm.core, title="")


## ----pca1, fig.cap="Primary Component Analysis of samples - all combinations of the 5 first primary components", out.width='100%'----
PCAplot.eset(affyNorm.core, categories1="genotype", categories2="treatment", title="", colorpalette=cbPalette, numcomponents=4)


## ----erccannots----------------------------------------------------------
annots <- read.csv(file.path(metaDir, "MoGene-2_0-st-v1.na34.mm10.transcript.csv"), skip=23)
erccindices <- which(grepl("ercc", annots$category))
erccprobesets <- annots$probeset_id[erccindices]


## ----loessnorm, results='hide'-------------------------------------------
eset.core <- exprs(affyNorm.core)
eset.corr.loess <- normalize.loess(eset.core, subset=erccindices, log.it=FALSE)
affyLoess.core <- affyNorm.core
exprs(affyLoess.core) <- eset.corr.loess


## ----cluster2, out.width='100%'------------------------------------------
plot_dendro(affyLoess.core, title="", labels.colname="treatment", colors.colname="genotype")


## ----PCAsd2, out.width='75%'---------------------------------------------
PCAplot.sd.eset(affyNorm.core, title="")


## ----pca2, fig.cap="Primary Component Analysis of samples - all combinations of the 5 first primary components", out.width='100%'----
PCAplot.eset(affyNorm.core, categories1="genotype", categories2="treatment", title="Genotypes", colorpalette=cbPalette, numcomponents=4)


## ----quantro-------------------------------------------------------------
library(doParallel)
cl <- makeCluster(16)
registerDoParallel(cl)

qtest <- quantro(object=exprs(affyRaw), groupFactor=as.vector(pData(affyRaw)$group), B=1000)
stopCluster(cl)
save.image(file.path(resultsDir, "RDATA"))

q()
