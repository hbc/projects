## ----setup---------------------------------------------------------------
if(file.exists("/Users/johnhutchinson/Work/projects")){
  baseDir  <- "/Users/johnhutchinson/Work/projects/aw_NMF_proteomics/June2015/"
  } else if (file.exists("/n/home08/")){
    baseDir <- "/n/hsphS10/hsphfs1/chb/projects/aw_NMF_proteomics/June2015"
    }
dataDir <- file.path(baseDir, "data")
resultsDir <- file.path(baseDir, "results")
metaDir <- file.path(baseDir, "meta")

## ----libraries-----------------------------------------------------------
library(reshape2)
library(ggplot2)
library(readr)
library(dplyr)
library(pheatmap)
library(ggdendro)
library(CHBUtils)
library(edgeR)
library(plyr)
library(RUVSeq)
library(limma)
library(knitr)
library(NMF)

## ----loaddata------------------------------------------------------------
rawdata.set0 <- read.csv(file.path(dataDir, "Set0", "Set0-Summary-2.csv"))
rawdata.set1 <- read.csv(file.path(dataDir, "Set1", "Set1-Summary-2.csv"))
rawdata.set2 <- read.csv(file.path(dataDir, "Set2", "Set2-Summary-2.csv"))

## ----summarizephospresidues----------------------------------------------
sumdata.set0 <- aggregate(cbind(Summed.126, Summed.127, Summed.128, Summed.129, Summed.130, Summed.131) ~ GeneName+Protein.Relative.Modifications.1, data=rawdata.set0, sum)
sumdata.set1 <- aggregate(cbind(Summed.126, Summed.127, Summed.128, Summed.129, Summed.130, Summed.131) ~ GeneName+Protein.Relative.Modifications.1, data=rawdata.set1, sum)
sumdata.set2 <- aggregate(cbind(Summed.126, Summed.127, Summed.128, Summed.129, Summed.130, Summed.131) ~ GeneName+Protein.Relative.Modifications.1, data=rawdata.set2, sum)

## ----dataclean-----------------------------------------------------------
# remove any rows with zero counts for all samples
sumdata.set0 <- sumdata.set0[!(apply(sumdata.set0[, grep("Summed", names(sumdata.set0))], 1, function(x) all(x==0))),]
sumdata.set1 <- sumdata.set1[!(apply(sumdata.set0[, grep("Summed", names(sumdata.set1))], 1, function(x) all(x==0))),]
sumdata.set2 <- sumdata.set2[!(apply(sumdata.set0[, grep("Summed", names(sumdata.set2))], 1, function(x) all(x==0))),]

## ----dataexplore---------------------------------------------------------
# merge datasets together - all phospho sites
sumdata.set.0.1. <- merge(sumdata.set0, sumdata.set1, by=c("GeneName", "Protein.Relative.Modifications.1"), suffixes=c(".0", ".1"), all=TRUE)
sumdata.set.0.1.2 <- merge(sumdata.set.0.1., sumdata.set2, by=c("GeneName", "Protein.Relative.Modifications.1"), all=TRUE)
sumdata <- sumdata.set.0.1.2
names(sumdata)[(ncol(sumdata)-5):ncol(sumdata)] <- paste(names(sumdata)[(ncol(sumdata)-5):ncol(sumdata)], ".2", sep="")
sumdata.m <- melt(sumdata, id.vars=c("GeneName", "Protein.Relative.Modifications.1"))
ggplot(sumdata.m, aes(x=value, col=variable) )+geom_density()+scale_x_log10()+ggtitle("Summed counts, all phospho sites")


# merge datasets together - only common phospho sites
sumdata.set.0.1. <- merge(sumdata.set0, sumdata.set1, by=c("GeneName", "Protein.Relative.Modifications.1"), suffixes=c(".0", ".1"))
sumdata.set.0.1.2 <- merge(sumdata.set.0.1., sumdata.set2, by=c("GeneName", "Protein.Relative.Modifications.1"))
sumdata <- sumdata.set.0.1.2
names(sumdata)[(ncol(sumdata)-5):ncol(sumdata)] <- paste(names(sumdata)[(ncol(sumdata)-5):ncol(sumdata)], ".2", sep="")
sumdata.m <- melt(sumdata, id.vars=c("GeneName", "Protein.Relative.Modifications.1"))
ggplot(sumdata.m, aes(x=value, col=variable) )+geom_density()+scale_x_log10()+ggtitle("Summed counts, common phospho sites")

## ----normvalues----------------------------------------------------------
# using just total count norm
## merge datasets together - only common phospho sites
sumdata.set.0.1. <- merge(sumdata.set0, sumdata.set1, by=c("GeneName", "Protein.Relative.Modifications.1"), suffixes=c(".0", ".1"))
sumdata.set.0.1.2 <- merge(sumdata.set.0.1., sumdata.set2, by=c("GeneName", "Protein.Relative.Modifications.1"))
sumdata <- sumdata.set.0.1.2
names(sumdata)[(ncol(sumdata)-5):ncol(sumdata)] <- paste(names(sumdata)[(ncol(sumdata)-5):ncol(sumdata)], ".2", sep="")
## calc total intensities for each sample
colsums <- colSums(sumdata[, grep("Summed", names(sumdata))])
## calc multiplier modifier for each samples, based on sample with total lowest intensity
mods<- 1/(colsums/min(colsums))
## normalize summed data with multiplier modifier
normed.sums <- as.data.frame(t(t(sumdata[,grep("Summed", names(sumdata))])*mods))
## add normalized data to non-normed data
names(normed.sums) <- sub("Summed", "Normed", names(normed.sums))
data <- cbind(sumdata, normed.sums)
normed.data <- cbind(data[,1:2], data[,grep("Normed", names(data))])
normed.data$gene_phosphosite <- paste(normed.data$GeneName, normed.data$Protein.Relative.Modifications.1, sep="_")
row.names(normed.data) <- normed.data$gene_phosphosite
normed.data <- log2(normed.data[,grep("Normed", names(normed.data))] + 0.5)

# plot normed results
normed.data.m <- melt(normed.data)

ggplot(normed.data.m, aes(x=value, col=variable) )+geom_density()
heatmap(as.matrix(normed.data), labRow = NA)
myDist <- dist(t(1-cor(normed.data)))
myTree <- hclust(myDist, method = "ward.D2")
dhc <- as.dendrogram(myTree)
ggdendrogram(dhc)

## ----edgeRnorm-----------------------------------------------------------
# merge datasets together - all phospho sites
sumdata.set.0.1. <- merge(sumdata.set0, sumdata.set1, by=c("GeneName", "Protein.Relative.Modifications.1"), suffixes=c(".0", ".1"), all=TRUE)
sumdata.set.0.1.2 <- merge(sumdata.set.0.1., sumdata.set2, by=c("GeneName", "Protein.Relative.Modifications.1"), all=TRUE)
sumdata <- sumdata.set.0.1.2
names(sumdata)[(ncol(sumdata)-5):ncol(sumdata)] <- paste(names(sumdata)[(ncol(sumdata)-5):ncol(sumdata)], ".2", sep="")

# prep for edgeR import
phosphosites <- paste(sumdata[,1], sumdata[,2], sep="_")
sumdata <- cbind(phosphosites, sumdata[, grep("umm", names(sumdata))])
row.names(sumdata) <- sumdata$phosphosites
sumdata$phosphosites <- NULL
sumdata[is.na(sumdata)] <- 0

# import into edgeR object and calc norm factors
sumdata.dge <- DGEList(counts=sumdata)
sumdata.dge <- calcNormFactors(sumdata.dge)

# output normed data adjusted for library size and
normed.data.edger <- cpm(sumdata.dge, normalized.lib.sizes = TRUE, log = TRUE)

# subset to commmon phospho sites
normed.data.edger <- normed.data.edger[row.names(normed.data.edger) %in% row.names(normed.data),]
normed.data.edger <- as.data.frame(normed.data.edger)

# plot normed results
normed.data.edger.m <- melt(normed.data.edger)

ggplot(normed.data.edger.m, aes(x=value, col=variable) )+geom_density()
heatmap(as.matrix(normed.data.edger), labRow = NA)
myDist <- dist(t(1-cor(normed.data.edger)))
myTree <- hclust(myDist, method = "ward.D2")
dhc <- as.dendrogram(myTree)
ggdendrogram(dhc)

## ----svaseq, eval=FALSE--------------------------------------------------
## metadata <- laply(strsplit(names(normed.data), "\\."), function(x) {
##   sample <- x[2]
##   batch <- x[3]
##   return(list(sample=sample, batch=batch))
## })
## metadata <- as.data.frame(metadata)
## metadata$sampleclass <- ifelse(metadata$sample==126| metadata$sample==127 | metadata$sample==128, "D1", "D4")
## 
## 
## set = newSeqExpressionSet(as.matrix(normed.data.edger))
## difference <- matrix(data=c(c(1:3,7:9,13:15), c(4:6,10:12,16:18)), byrow=TRUE, nrow=2)
## batch_ruv_emp <- RUVs(as.matrix(normed.data.edger), rownames(normed.data.edger), k=2, difference)
## 
## normed.suv <- as.data.frame(log2(batch_ruv_emp$normalizedCounts + 0.5))
## 
## normed.suv.m <- melt(normed.suv)
## 
## ggplot(normed.suv.m, aes(x=value, col=variable) )+geom_density()
## heatmap(as.matrix(normed.suv), labRow = NA)
## myDist <- dist(t(1-cor(normed.suv)))
## myTree <- hclust(myDist, method = "ward.D2")
## dhc <- as.dendrogram(myTree)
## ggdendrogram(dhc)

## ----loadbatchcorrected--------------------------------------------------
load(file.path(resultsDir, "batch.corrected.counts"))

## ----limma, results='asis'-----------------------------------------------
row.names(metadata) = names(sumdata.dge$counts)
dd = cbind(metadata, batch_ruv_emp$W)[,1:5]
dd$batch = as.factor(unlist(dd$batch))  

ma = model.matrix(~ 0 + sampleclass + batch, data=dd)
  
dat.voom = voom(batch_ruv_emp$normalizedCounts,design = ma, plot = TRUE)
dat.fit <- lmFit(dat.voom, ma)

cont.ma = makeContrasts(condition=sampleclassD4-sampleclassD1, levels=ma)
dat.fit.cont <- contrasts.fit(dat.fit, cont.ma)
dat.bayes <- eBayes(dat.fit.cont)

kable(summary(decideTests(dat.bayes)))

all_de = topTable(dat.bayes, coef="condition", number = Inf)

## ----output--------------------------------------------------------------
write.csv(normed.suv, file=file.path("normalized.data.suvseq.csv"))
write.csv(all_de, file=file.path("limma.voom.csv"))
# write.csv(normed.data.edger, file=file.path(dataDir, "Set2", "normalized.data.edgeR.csv"))

## ----nmfvars-------------------------------------------------------------
TRAIN=50
RUN=1250
mad.cutoff=0.5

minnumfeatures=25

## ----estimrank-----------------------------------------------------------
# drop weird rows with all negative numbers or only zeroes
eset.corr <- normed.suv
eset.corr <- eset.corr[apply(eset.corr, 1, function(x) all(x>0)),]

groups.corr <-  as.factor(metadata$sampleclass)

estim.corr <- nmf(eset.corr, 2:5, nrun = TRAIN, seed = 123456, .options='v-p') #disable parallel compute
plot(estim.corr)

## ----overfitcheck, results='hide',warning=FALSE, message=FALSE, error=FALSE,cache=TRUE----
# shuffle original data to look for overfitting
eset.corr.rand <- randomize(eset.corr)
# estimate quality measures from the shuffled data (use default NMF algorithm)
estim.corr.rand <- nmf(eset.corr.rand, 2:5, nrun = TRAIN, seed = 12345, .options="v-p")
# plot measures on same graph
plot(estim.corr, estim.corr.rand)

## ----estimatefactoriziationrank.qualitative, results='hide'--------------
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

consensusmap(estim.corr, annCol=list(samplegroups=as.character(unlist(eset.corr$group_short ))),annColors=list(samplegroups=cbPalette[1:3]), labCol=groups.corr, labRow=groups.corr, scale="row", color='-RdYlBu2:200')

## ----comparealgs, cache=TRUE, results="hide", message=FALSE, error=FALSE,warning=FALSE----
res.multi.method.2 <- nmf(eset.corr, 2, list("brunet", "KL", "lee","nsNMF"), nrun=TRAIN, seed = 123456, .options = "tv-p")
plot(res.multi.method.2, main="NMF residuals - 2 metagenes")

## ----fullNMF, eval=TRUE--------------------------------------------------
res.final.2 <- nmf(eset.corr, 2, "nsNMF", nrun=RUN, .options = "tv-p")

# save precomputed NMF values, hack to avoid using caching
save(list="res.final.2", file=file.path(resultsDir,  "RDATA.res.final.2"))

