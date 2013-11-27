## @knitr libraries1
library(affy)
library(simpleaffy)
library(arrayQualityMetrics)
library(limma)
library(mouse430a2.db)
library(pheatmap)
library(RColorBrewer)
library(xtable)
library(biomaRt)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
library(plyr)
library(ggplot2)
library(reshape)
library(methylKit)
library(parallel)
library(IRanges)
library(GenomicRanges)
library(gpairs)


## @knitr functions1
rbind.fill <- function(l) {
  do.call(rbind, lapply(lapply(l, unlist), "[", unique(unlist(c(sapply(l,names))))))
  }

bedTools.2in<-function(functionstring="intersectBed",bed1,bed2,opt.string="")
{
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) ## do not not to use scientific notation when writing out
   ## write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
   ## create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  print("pasg")
  cat(command,"\n")
  try(system(command))
   res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}

getCorr = function(obj){
   meth.mat  <-  getData(obj)[,obj@numCs.index]/(obj[,obj@numCs.index]+obj[,obj@numTs.index])  
   names(meth.mat) <- obj@sample.ids
   return(meth.mat)
}


## @knitr datadirectories
baseDir="/n/hsphS10/hsphfs1/chb/projects/JH/kc_rrbs"
annotDir="/n/hsphS10/hsphfs1/chb/projects/JH/kc_rrbs/data/annotations"
exprdataDir=paste(baseDir, "data/RNA_Microarray", sep="/")
methdataDir <- paste(baseDir, "data/fastq/trimmed", sep="/")
resultsDir=paste(baseDir, "results", sep="/")


## @knitr exprdataload
mic.data.raw <- read.affy('covars.desc', path=exprdataDir, verbose=T)


## @knitr metadata
pData(mic.data.raw)$group <- paste(pData(mic.data.raw)$origin , pData(mic.data.raw)$tissue, sep="_")
pd <- pData(mic.data.raw)
metadataTable <- xtable(pd)
print(metadataTable, type='html')
rm(metadataTable) ## cleanup


## @knitr normalize_and_extract_expression_values
mic.data.norm <- call.exprs(mic.data.raw, "rma")
eset.norm <- exprs(mic.data.norm)
colnames(eset.norm) <- pData(mic.data.norm)$Sample


## @knitr rawQC
arrayQualityMetrics(expressionset=mic.data.raw, outdir=paste(resultsDir, "report_raw", sep="/"), force=TRUE, do.logtransform=TRUE)
rm(mic.data.raw) ## cleanup


## @knitr normQC
arrayQualityMetrics(expressionset=mic.data.norm, outdir=paste(resultsDir, "report_rma", sep="/"), force=TRUE, do.logtransform=FALSE)


## @knitr exprPCA
myPca <- prcomp(t(eset.norm))
# SD of components
plot(myPca$sdev, xlab="PCA", ylab="sddev", main="Variance explained by each principal component")
# Plot samples in all combinations in the 5 first primary components
tmpPCAData <- as.data.frame(myPca$x[,1:5])
colors <- brewer.pal(length(unique(paste(pData(mic.data.norm)$tissue, pData(mic.data.norm)$origin, sep="_"))), "Set1")[factor(paste(pData(mic.data.norm)$tissue, pData(mic.data.norm)$origin, sep="_"))]
plot(tmpPCAData, col=colors, pch=row.names(myPca$x), main="PCA plot of first 5 components")
rm(myPca, tmpPCAData, colors) ## cleanup


## @knitr annotations
symbols=mget(row.names(eset.norm), mouse430a2SYMBOL, ifnotfound=NA)
ensembls=mget(row.names(eset.norm),mouse430a2ENSEMBL, ifnotfound=NA)
gene.annots=as.data.frame(cbind(symbols, ensembls))
nrow(gene.annots)==nrow(eset.norm) ## quality check, if TRUE, annotations are good
all(row.names(eset.norm)==row.names(gene.annots)) ## another quality check, if TRUE, annotations are good


## @knitr mediansummarize
eset.genes.norm <- aggregate(eset.norm, by=list(symbols=as.vector(unlist(gene.annots$symbols))), median)
row.names(eset.genes.norm) <- eset.genes.norm[,"symbols"]
eset.genes.norm <- eset.genes.norm[,-(grep("symbols", names(eset.genes.norm)))]


## @knitr design
design <- model.matrix(~ -1+factor(pData(mic.data.norm)$group))
## always make sure the headings match
colnames(design) <- c("ES_ES", "iPS_TTF", "iPS_VM", "somatic_TTF", "somatic_VM")

designTable <- xtable(design)
print(designTable, type='html')
rm(designTable) ## cleanup


## @knitr contrastmatrix
design.pairs <- function(levels) {
 n <- length(levels)
 design <- matrix(0, n, choose(n, 2))
 rownames(design) <- levels
 colnames(design) <- 1:choose(n, 2)
 k <- 0
 for (i in 1:(n - 1))
   for (j in (i + 1):n) {
     k <- k+1
     design[i, k] <- 1
     design[j, k] <- -1
     colnames(design)[k] <- paste(levels[i], "-", levels[j], sep="")
     }
 design
}
contrast.matrix <- design.pairs(levels(as.factor(pData(mic.data.norm)$group)))

contrastmatrixTable <- xtable(contrast.matrix)
print(contrastmatrixTable, type='html')
rm(contrastmatrixTable, design.pairs) ## cleanup


## @knitr linearmodel
fit <- lmFit(eset.genes.norm, design) 


## @knitr contrastfit
fit2 <- contrasts.fit(fit, contrast.matrix)



## @knitr bayes
fit2 <- eBayes(fit2) 


## @knitr balance
diffs.temp=as.data.frame(rbind.fill(apply(decideTests(fit2,method="separate",  lfc=4),2,table)))
names(diffs.temp)=c("0","1","-1")
diffs.temp$log2ratio=log2(diffs.temp[,"1"]/diffs.temp[,"-1"])
balanceTable <- xtable(diffs.temp)
print(balanceTable, type='html')
rm(balanceTable, diffs.temp) ## cleanup


## @knitr volcanoplot
llply(setdiff(seq(1,ncol(fit2$contrasts), 1),grep("somatic", dimnames(fit2$contrasts)[[2]])), function(n){
  n <- as.numeric(unlist(n))
  stats.core <- topTable(fit2, coef=n, sort.by="B",number=length(symbols), genelist=fit2$genes)
  stats.core$Passes.0.05_FDR.threshold  <-  as.factor(stats.core$adj.P.Val<0.05)
  g <- ggplot(data=stats.core, aes(x=logFC, y=-log10(P.Value), color=Passes.0.05_FDR.threshold, size=B)) +
    geom_point(alpha=0.5) + 
    geom_vline(xintercept=c(-1,1), color="orange", alpha=0.8) +
    xlab("log2 fold change") + ylab("-log10 p-value") +
    labs(title=dimnames(fit2$contrasts)[[2]][n])
  return(g)
  })


## @knitr allexprstats
all.stats=lapply(seq(1,ncol(contrast.matrix),1), function(x) {
  contrast=colnames(contrast.matrix)[x]
  stats.all.core=topTable(fit2, coef=x, adjust="fdr", p.value=1,sort="B",number=nrow(fit2$genes), genelist=fit2$genes)
  eset.all.core  <-  eset.genes.norm[stats.all.core[, "ID"], ]
  stats=cbind(stats.all.core, eset.all.core)
  return(list(stats=stats, contrast=contrast))
})
names(all.stats)=dimnames(fit2$contrasts)[[2]]


## @knitr outputs
lapply(all.stats, function(x) {
  contrast=x$contrast
  write.table(x$stats, file=paste("../results/all.genes.stats", contrast, "xls",sep="."), sep="\t", row.names=F, col.names=T)
  })


## @knitr methsampleinfo
sampleinfo <- read.table(paste(methdataDir, "sampleIDs.tab", sep="/"), header=T)
datafiles <- list.files(methdataDir, pattern="CpG.txt")
sampleinfo$files <- datafiles[as.vector(apply(sapply(as.vector(sampleinfo$fileprefix), regexpr, datafiles), 2, function(n) which(n==1)))]


## @knitr methdataimport
sampleIDs <- as.list(as.vector(sampleinfo$sampleID))
locations <- as.list(paste(methdataDir, sampleinfo$files, sep="/"))
treatments <- as.vector(sampleinfo$treatment)
meth.quants <-read(location=locations, sample.id=sampleIDs, assembly="mm9", context="CpG", pipeline="bismark", resolution="base", treatment=treatments) ##long step
rm(sampleIDs, locations, treatments, methdataDir)


## @knitr methstats1
lapply(meth.quants, function(n) {
  print(n@sample.id)
  getMethylationStats(n, plot=FALSE)
  })
lapply(meth.quants, function(n) getMethylationStats(n, plot=TRUE, labels=FALSE))


## @knitr methstats2
lapply(meth.quants, function(n) getCoverageStats(n, plot=TRUE, labels=TRUE))


## @knitr methcoveragefilter
meth.quants.filtered <- filterByCoverage(meth.quants, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
rm(meth.quants) ## cleanup


## @knitr methmerge
meth.quants.filtered.merged <- unite(meth.quants.filtered, destrand = TRUE)
meth.quants.filtered.merged@sample.ids=sub("_p8", "", sub("_derived", "", sub("Mouse_", "", meth.quants.filtered.merged@sample.ids)))


## @knitr correlations
cor.mat=getCorr(meth.quants.filtered.merged)
gpairs(cor.mat[sample(1:nrow(cor.mat),10000 ),], lower.pars=list(scatter="lm"), upper.pars=list(scatter="stats"), stat.pars=list(fontsize=16, use.color=FALSE), scatter.pars=list(pch=20, col="#00000022"), diag.pars=list(fontsize=10, show.hist=FALSE))
rm(cor.mat) ## cleanup


## @knitr clustering
meth.dendrogram=clusterSamples(meth.quants.filtered.merged, dist="correlation", method="ward", plot=FALSE)
plot(meth.dendrogram)
rm(meth.dendrogram) ## cleanup


## @knitr PCA.screeplot
PCASamples(meth.quants.filtered.merged, screeplot = TRUE)


## @knitr PCA.dimension.plot
PCASamples(meth.quants.filtered.merged)


## @knitr tile
tiles <- tileMethylCounts(meth.quants.filtered, win.size = 1000, step.size = 1000)
tiles.merged <- unite(tiles)
tiles.merged@sample.ids=sub("_p8", "", sub("_derived", "", sub("Mouse_", "", tiles.merged@sample.ids)))
rm(tiles) ## cleanup


## @knitr functionalregions
CpGIslands.bed=read.delim(paste(annotDir, "mm9_CpG_Islands.bed", sep="/"), header=F)
names(CpGIslands.bed)=c("chr", "start", "end", "id")
## subset to chromosomes present in methylation reads
CpGIslands.bed=CpGIslands.bed[CpGIslands.bed$chr %in% meth.quants.filtered[[1]]@.Data[[2]],]
CpGIslands.bed=droplevels(CpGIslands.bed)
## convert to GRanges objects
CpGIslands.GR=GRanges(seqnames=CpGIslands.bed$chr, ranges=IRanges(start=CpGIslands.bed$start, end=CpGIslands.bed$end), ids=CpGIslands.bed$id)
## get methylation counts for regions
meth.quants.CpGIslands=regionCounts(meth.quants.filtered,  CpGIslands.GR)
rm(meth.quants.filtered, CpGIslands.GR) ## cleanup
## merge 
meth.quants.CpGIslands.merged <- unite(meth.quants.CpGIslands)
## change ids
meth.quants.CpGIslands.merged@sample.ids=sub("_p8", "", sub("_derived", "", sub("Mouse_", "", meth.quants.CpGIslands.merged@sample.ids)))



## @knitr plot_percent_meth_dists
data.perc.meth=meth.quants.CpGIslands.merged[,c("id", "chr", "start", "end","strand")]
for (n in 1:length(grep("numCs", names(meth.quants.CpGIslands.merged)))) {
  coverage=meth.quants.CpGIslands.merged[,paste("coverage",n, sep="")]
  numCs=meth.quants.CpGIslands.merged[,paste("numCs",n, sep="")]
  perc.meth=numCs/coverage
  data.perc.meth=cbind(data.perc.meth, perc.meth) 
  names(data.perc.meth)[ncol(data.perc.meth)]=sub("$", n,  names(data.perc.meth)[ncol(data.perc.meth)])
} 
rm(coverage, numCs, perc.meth)
data.perc.meth.m=melt(data.perc.meth, id.vars=c("id", "chr", "start", "end", "strand"))

ggplot(data.perc.meth.m, aes(x=value, color=variable))+geom_density() +
  scale_x_log10() +
  xlab("% methylation")+
  scale_color_hue( name="Sample",  breaks=c("perc.meth1","perc.meth2","perc.meth3","perc.meth4","perc.meth5","perc.meth6", "perc.meth7", "perc.meth8", "perc.meth9"),labels=unlist(lapply(meth.quants.CpGIslands, function(n) n@sample.id)))



## @knitr diffmeth.tiled
## tiled regions
diffs.all.tiled=alply(combos, 2, function(n) {
  indices=which(tiles.merged@treatment %in% n)
  sample.ids.subset=tiles.merged@sample.ids[indices]
  treatments.subset=c(0,1)[factor(tiles.merged@treatment[indices])]
  tiles.merged.subset=reorganize(tiles.merged, sample.ids=sample.ids.subset,treatment=treatments.subset)
  diffs=calculateDiffMeth(tiles.merged.subset, num.cores=detectCores())
  return(diffs)
  })
names(diffs.all.tiled)=combos.by.name


## @knitr diffmeth.CpGs
diffs.all.CpGIslands=alply(combos, 2, function(n) {
  indices=which(meth.quants.CpGIslands.merged@treatment %in% n)
  sample.ids.subset=meth.quants.CpGIslands.merged@sample.ids[indices]
  treatments.subset=c(0,1)[factor(meth.quants.CpGIslands.merged@treatment[indices])]
  meth.quants.CpGIslands.subset=reorganize(meth.quants.CpGIslands.merged, sample.ids=sample.ids.subset, treatment=treatments.subset)
  diffs=calculateDiffMeth(meth.quants.CpGIslands.subset, num.cores=detectCores())
  return(diffs)
  })
names(diffs.all.CpGIslands)=combos.by.name


## @knitr diffmethresults.tiled
lapply(diffs.all.tiled, function(n) {
  diffs.all.tiled[[1]]@sample.id
diffMethPerChr(n, meth.cutoff=25, qvalue=0.05, plot=TRUE)
})
myDiffs25p.tiled <- lapply(diffs.all.tiled, function(x) get.methylDiff(x, difference = 25,    qvalue = 0.01))
myDiffs25p.hypo.tiled <- lapply(diffs.all.tiled, function(x) get.methylDiff(x, difference = 25,    qvalue = 0.01, type="hypo"))
myDiffs25p.hyper.tiled <- lapply(diffs.all.tiled, function(x) get.methylDiff(x, difference = 25,    qvalue = 0.01, type="hyper"))


## @knitr diffmethresults.CpGIslands
lapply(diffs.all.CpGIslands, function(n) {
  print(diffs.all.CpGIslands[[1]]@sample.id)
diffMethPerChr(n, meth.cutoff=25, qvalue=0.05, plot=TRUE)
})
myDiffs25p.CpGI <- lapply(diffs.all.CpGIslands, function(x) get.methylDiff(x, difference = 25,    qvalue = 0.01))
myDiffs25p.hypo.CpGI <- lapply(diffs.all.CpGIslands, function(x) get.methylDiff(x, difference = 25,    qvalue = 0.01, type="hypo"))
myDiffs25p.hyper.CpGI <- lapply(diffs.all.CpGIslands, function(x) get.methylDiff(x, difference = 25,    qvalue = 0.01, type="hyper"))



## @knitr get_genomic_features
gene.obj=read.transcript.features(paste(exprdataDir, "/mm9.knowngenes.bed", sep=""))
CpGI.obj=read.feature.flank(paste(exprdataDir, "/mm9CpGIslands.bed", sep=""), feature.flank.name = c("CpGi","shores"))
lapply(myDiffs25p.tiled, function(n) annotate.WithGenicParts(n, gene.obj))
lapply(myDiffs25p.tiled, function(n) annotate.WithFeature.Flank(n, CpGI.obj$CpGi, CpGI.obj$shores, feature.name = "CpGi",flank.name = "shores"))


## @knitr assign_CpGIs_to_genes
known_genes.bed=read.delim(paste(annotDir, "mm9_known_genes.bed", sep="/"), header=F)[,c(1:6)]
promoters.bed=read.delim(paste(annotDir,"5kb_upstream_of_mm9_known_genes.bed", sep="/"), header=F)[,c(1:6)]
names(promoters.bed)=c("chr", "start", "end", "id", "score", "strand")
names(known_genes.bed)=c("chr", "start", "end", "id", "score", "strand")
## convert ensembl transcriptids to ensembl geneids and remake bed file
known_genes.bed=merge(getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id"), filters="ensembl_transcript_id", value=known_genes.bed$id, mart=ensembl), known_genes.bed, by.x="ensembl_transcript_id", by.y="id")
known_genes.bed=known_genes.bed[,c(3,4,5,2,6,7)]
promoters.bed=merge(getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id"), filters="ensembl_transcript_id", value=promoters.bed$id, mart=ensembl), promoters.bed, by.x="ensembl_transcript_id", by.y="id")
promoters.bed=promoters.bed[,c(3,4,5,2,6,7)]
## relate CpG islands to promoters
CpGIslands.promoters<-bedTools.2in( bed1=CpGIslands.bed, bed2=promoters.bed, opt.string="-wo -f 1")
names(CpGIslands.promoters)<-c("chr.CpGI", "start.CpGI", "end.CpGI", "id.CpGI", "chr.prom", "start.prom", "end.prom", "ensembl_gene_id.prom", "score.prom", "strand.prom", "bp.overlap")
CpGIslands.promoters$CpGI.loc.ID=paste(CpGIslands.promoters$chr.CpGI, CpGIslands.promoters$start.CpGI, CpGIslands.promoters$end.CpGI, sep=".")


## @knitr plot_relationships
## merge CpG Island methylation and gene expression data for all comparisons through "CpGIslands.promoters" table
lapply(as.list(names(all.stats)), function(n) {
  stats.exprs=all.stats[[n]]$stats
  stats.exprs=merge(stats.exprs, gene.annots, by.x="ID", by.y="symbols", all=FALSE)
  stats.CpGmeth=getData(diffs.all.CpGIslands[[n]])
  names(stats.CpGmeth)=paste(names(stats.CpGmeth), "CpGI", sep=".")
  stats.CpGmeth=merge(stats.CpGmeth, CpGIslands.promoters, by.x=c("id.CpGI","chr.CpGI", "start.CpGI", "end.CpGI"),  by.y=c("CpGI.loc.ID", "chr.CpGI", "start.CpGI", "end.CpGI"), all=FALSE)
  stats=merge(stats.CpGmeth, stats.exprs, by.x="ensembl_gene_id.prom", by.y="ensembls", all=FALSE)
ggplot(data=subset(stats, qvalue.CpGI<0.05 & adj.P.Val<0.05), aes(x=logFC, y=meth.diff.CpGI, color=-log10(P.Value), size=-log10(pvalue.CpGI), label=ID))+geom_point(alpha=0.4)+scale_color_gradient(low="blue", high="red")+scale_size_continuous(range=c(2,12))+geom_text(size=5, vjust=-1.1)+geom_hline(yintercept=0, color="darkgrey")+geom_vline(xintercept=0, color="darkgrey")+labs(title=n)
  })


## @knitr relate_somatic_differences_to_embryonic
#bring together expression data for all comparisons into one large table
all.merged.stats.exprs=all.stats[[1]]$stats[,1:7]
for (n in 2:length(names(all.stats))){
  all.merged.stats.exprs=merge(all.merged.stats.exprs, all.stats[[n]]$stats[,1:7], by="ID", suffixes=c(paste(".", names(all.stats)[n-1], sep=""), paste(".", names(all.stats)[n], sep="")))
  }
all.merged.stats.exprs=merge(all.merged.stats.exprs, gene.annots, by.x="ID", by.y="symbols", all=FALSE)
#bring together methylation data for all comparisons into one large table
all.merged.stats.CpGI=diffs.all.CpGIslands[[1]]
for (n in 2:length(names(all.stats))){
  all.merged.stats.CpGI=merge(all.merged.stats.CpGI, diffs.all.CpGIslands[[n]], by=c("id", "chr","start","end","strand"), suffixes=c(paste(".", names(all.stats)[n-1], sep=""), paste(".", names(all.stats)[n], sep="")))
  }
names(all.merged.stats.CpGI)=sub("$", ".CpGI", names(all.merged.stats.CpGI))
all.merged.stats.CpGI=merge(all.merged.stats.CpGI, CpGIslands.promoters, by.x=c("id.CpGI","chr.CpGI", "start.CpGI", "end.CpGI"),  by.y=c("CpGI.loc.ID", "chr.CpGI", "start.CpGI", "end.CpGI"), all=FALSE)
#merge expression and methylation data for all comparisons into one large table
all.merged.stats=merge(all.merged.stats.CpGI, all.merged.stats.exprs, by.x="ensembl_gene_id.prom", by.y="ensembls", all=FALSE)
names(all.merged.stats)=sub("-", "_", names(all.merged.stats))



#plot relationship between somatic methylation differences and iPS methylation differences
ggplot(subset(all.merged.stats, qvalue.iPS_TTF_iPS_VM.CpGI<0.01 & qvalue.somatic_TTF_somatic_VM.CpGI<0.01), aes(x=meth.diff.iPS_TTF_iPS_VM.CpGI, y=meth.diff.somatic_TTF_somatic_VM.CpGI, size=-log10(qvalue.iPS_TTF_iPS_VM.CpGI), color=-log10(qvalue.somatic_TTF_somatic_VM.CpGI), label=ID))+
  geom_point()+ 
  scale_color_gradient(low="blue", high="red", guide="colourbar")+
  scale_size_continuous(range=c(2,12))+
  geom_text(size=5, vjust=-1.1)+
  ggtitle("CpGIsland Methylation differences \n between TTF and VM samples in both somatic and iPS cells")+
  ylab("somatic")+xlab("iPS")+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)


#plot relationship between ES-iPS_TTF methylation differences and ES-somatic_TTF methylation differences
ggplot(subset(all.merged.stats, qvalue.ES_ES_iPS_TTF.CpGI<0.01 & qvalue.ES_ES_iPS_VM.CpGI<0.01), aes(x=meth.diff.ES_ES_iPS_TTF.CpGI, y=meth.diff.ES_ES_iPS_VM.CpGI, size=-log10(qvalue.ES_ES_iPS_TTF.CpGI), color=-log10(qvalue.ES_ES_iPS_VM.CpGI), label=ID))+
  geom_point()+ 
  scale_color_gradient(low="blue", high="red", guide="colourbar")+
  scale_size_continuous(range=c(2,12))+
  geom_text(size=5, vjust=-1.1)+
  ggtitle("CpGIsland Methylation differences \n between iPS/somatic TTF and ES cells")+
  ylab("ES vs. iPS_VM")+xlab("ES vs. iPS_TTF")+  
  geom_hline(yintercept=0)+geom_vline(xintercept=0)

#plot relationship between ES-iPS_VM methylation differences and ES-somatic_VM methylation differences
ggplot(subset(all.merged.stats, qvalue.ES_ES_iPS_TTF.CpGI<0.01 & qvalue.ES_ES_somatic_TTF.CpGI<0.01), aes(x=meth.diff.ES_ES_iPS_TTF.CpGI, y=meth.diff.ES_ES_somatic_TTF.CpGI, size=-log10(qvalue.ES_ES_iPS_TTF.CpGI), color=-log10(qvalue.ES_ES_somatic_TTF.CpGI), label=ID))+
  geom_point()+ 
  scale_color_gradient(low="blue", high="red", guide="colourbar")+
  scale_size_continuous(range=c(2,12))+
  geom_text(size=5, vjust=-1.1)+
  ggtitle("CpGIsland Methylation differences \n between iPS/somatic VM and ES cells")+
  ylab("ES vs. somatic_TTF")+xlab("ES vs. iPS_TTF")+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)
save.image()


## @knitr sessioninfo
print(sessionInfo())
}


