
## setup report details
clientname="Serge Gregoire"
clientemail="SGREGOIRE2@PARTNERS.ORG"
labPI="Wu"
lablocation="MGH/HSCI"
analystname="John Hutchinson"
analystemail="jhutchin@hsph.harvard.edu"


## @knitr libraries, echo=TRUE
library(oligo)
library(limma)
library(xtable)
library(Biobase)
library(pd.mogene.1.0.st.v1)
library("mogene10sttranscriptcluster.db")
library(plyr)
library(ggplot2)
library(biomaRt)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
filters=listFilters(ensembl)
attributes=listAttributes(ensembl)
library(Rcade)
library(parallel)


## @knitr VARIABLES, echo=TRUE
## Setup Data and Results directory variables
if(file.exists("/n/home08/jhutchin/")){
  baseDir="/n/hsphS10/hsphfs1/chb/projects/sw_cardiomyocyte_differentiation"
  } else if (file.exists("/Volumes/ody/")){
  baseDir="/Volumes/ody/projects/sw_cardiomyocyte_differentiation"
}
metaDir=file.path(baseDir, "meta")
mic.dataDir=file.path(baseDir, "data/microarray")
cs.bamDir=file.path(baseDir, "results/chipseq/bowtie")
mic.resultsDir=file.path(baseDir, "results/microarray")
cs.resultsDir=file.path(baseDir, "results/chipseq")
mic.covarsfilename="covars.desc" # do not use full path
cs.covarsfilename="chipseq.covars.desc"

mic.grouplabel="treatment"
mic.samplelabel="sample"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


## @knitr dataload, results='hide', cache=TRUE
covars <- read.delim(file.path(metaDir, mic.covarsfilename), row.names=1) ## simple tab delimited file with CEL file in first column (no heading for this column) and sample metadata (i.e. sampleID, treatment group, batch etc.) in subsequent columns
celFiles <- list.celfiles(mic.dataDir, full.names=TRUE)
affyRaw <- read.celfiles(celFiles, pkgname="pd.mogene.1.0.st.v1")
pData(affyRaw) <- covars 
validObject(affyRaw) ## sanity check


## @knitr covars, results='asis'
## Sample information table
pDataTable <- xtable(pData(affyRaw))
print(pDataTable, type='html')


## @knitr normalize, results='hide'
affyNorm.core <- rma(affyRaw, target="core", background=TRUE, normalize=TRUE)
rm(affyRaw) # cleanup


## @knitr features, results='hide'
# retrieve NetAffx Biological Annotation
featureData(affyNorm.core) <- getNetAffx(affyNorm.core, "transcript")
symbols <-  unlist(mget(as.character(pData(featureData(affyNorm.core))$transcriptclusterid), mogene10sttranscriptclusterSYMBOL, ifnotfound=NA))
entrezids <- unlist(mget(as.character(pData(featureData(affyNorm.core))$transcriptclusterid), mogene10sttranscriptclusterENTREZID, ifnotfound=NA))

# check to make sure data is correct
identical(length(featureData(affyNorm.core)$probesetid), length(symbols)) # sanity check, sane=TRUE
identical(length(featureData(affyNorm.core)$probesetid), length(entrezids)) # sanity check, sane=TRUE
gene.annots <- as.data.frame(cbind(symbols, entrezids))
head(gene.annots$symbols[!is.na(gene.annots$symbols)]) # sanity check, sane=>see gene ids


## @knitr design, results="asis"
design <- model.matrix(~ -1+factor(pData(affyNorm.core)[,mic.grouplabel]))
# make sure the headings match
colnames(design) <- sub("factor.pData.affyNorm.core... mic.grouplabel..", "", colnames(design))
designTable <- xtable(design)
print(designTable, type='html')


## @knitr contrastmatrix, results='asis'
contrast.matrix <- makeContrasts(control_GFP-minus_dox,control_GFP-plus_dox, plus_dox-minus_dox, levels=c("control_GFP", "minus_dox", "plus_dox"))
contrastmatrixTable <- xtable(contrast.matrix)
print(contrastmatrixTable, type='html')


## @knitr linearmodel
eset.core <- exprs(affyNorm.core) 
fit.core <- lmFit(eset.core, design) 


## @knitr contrastfit
fit2.core <- contrasts.fit(fit.core, contrast.matrix) 


## @knitr bayes
fit2.core <- eBayes(fit2.core) 


## @knitr allstats
all.stats <- llply(seq(1,3,1), function(n) {
    contrast <- gsub(" ", "", dimnames(fit2.core$contrasts)$Contrasts[n])
    stats.core <- topTable(fit2.core, coef=n, sort.by="B",number=length(symbols), genelist=cbind(gene.annots[,c("symbols", "entrezids")], fit2.core$genes))
    return(list(stats.core=stats.core, contrast=contrast))
    })



## @knitr load_chipseq
cs.targets <- read.delim(file.path(metaDir, cs.covarsfilename), as.is = TRUE)
#cs.targets$filepath <- file.path(baseDir, cs.targets$filepath)


## @knitr gene_annotations
tss.anno <- getBM(attributes= c("entrezgene", "chromosome_name","transcript_start", "transcript_end", "strand"), mart= ensembl)
tss.anno <- tss.anno[order(tss.anno$chromosome_name),]
colnames(tss.anno) <- c("ENTREZ","chr","start","end","str")


## @knitr tss_bins
ChIPannoZones <- defineBins(tss.anno, zone=c(-1500, 1500), geneID="ENTREZ")


## @knitr Rcade
cl <- makeCluster(2, "SOCK")
#all.stats[[3]]$contrast
DE <- all.stats[[3]]$stats.core
DElookup <- list(GeneID="entrezids", logFC="logFC", B="B", "symbols" )
Rcade <- RcadeAnalysis(DE, ChIPannoZones, ChIPtargets=cs.targets, ChIPfileDir = cs.bamDir, shift = 0, cl=cl, DElookup=DElookup)
setwd(baseDir)

save.image("Feb5.Rdata")
q()


