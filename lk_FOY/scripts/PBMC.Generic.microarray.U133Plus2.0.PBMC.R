## ----general_libraries---------------------------------------------------
source("http://bioconductor.org/biocLite.R") # BioConductor script necessary for installing new BioC libraries with biocLite()
require(xtable) # table generation for reports
require(plyr) # library for iteratively working with data structures
require(ggplot2) # library for plotting 
require(RColorBrewer) # library for generating color palettes
require(googleVis) # library for presenting tables
source("http://dl.dropboxusercontent.com/u/4253254/Resources/functions.r")


## ----general_directories-------------------------------------------------
if (file.exists("/n/hsphS10/hsphfs1/chb/projects/lk_FOY")) {
  baseDir <- "/n/hsphS10/hsphfs1/chb/projects/lk_FOY"
  }  else if (file.exists("/Volumes/home08/jhutchin/consults/lk_FOY/")) {
    baseDir <- "/Volumes/home08/jhutchin/consults/lk_FOY"
    } else {
      baseDir <- "/Volumes/ody/consults/lk_FOY"
      }
dataDir <- file.path(baseDir, "data")
resultsDir <- file.path(baseDir, "results", "PBMC", "U133Plus2")
metaDir <- file.path(baseDir, "meta", "PBMC")


## ----subset_data---------------------------------------------------------
refined.metadata <- read.delim(file.path(metaDir,"unified.metadata.refined.PBMC.tab"))
U133.Plus2.data <- refined.metadata[which(!is.na(refined.metadata$age) & !is.na(refined.metadata$gender) & !is.na(refined.metadata$CEL_regex) & grepl("GPL570|A-AFFY-44", refined.metadata$platform)),]
write.table(U133.Plus2.data, file.path(metaDir, "unified.metadata.refined.PBMC.U133Plus2.0.tab"), quote=F, sep="\t", row.names=F, col.names=T)
write.table(U133.Plus2.data, file=file.path(metaDir, "unified.metadata.refined.PBMC.U133Plus2.0.xls"), sep="\t", row.names=F, col.names=T)


## ----microarray_analysis_libraries---------------------------------------
# to parse the CEL files and work with intensity values
library(affy) 
# for QC reports
library(arrayQualityMetrics)
# library to do stats 
library(limma) 
# pretty heatmaps
library(pheatmap) 
# annotations for the hgU1332.0Plus array
library(hgu133plus2.db) 


## ----microarray_analysis_variables---------------------------------------
# colorblind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
# age ranges to compare
child.age.range <- c(5,12)
adult.age.range <- c(18,40)
neonate.age.range <- c(0,4)
# THESE CANNOT OVERLAP!


## ----covariatedataframe--------------------------------------------------
# U1332.0Plus only for now
covartemplate.file <- "unified.metadata.refined.PBMC.U133Plus2.0.tab" 
covartemplate <- read.table(file.path(metaDir,covartemplate.file), header=T, colClasses="character", sep="\t")
# convert character ages to numeric to allow numeric subsetting 
covartemplate$age <- as.numeric(covartemplate$age) 
## limit samples to those within the described child, adult and neonate age ranges
covartemplate <- covartemplate[union(union(which(covartemplate$age<=max(child.age.range) & covartemplate$age>=min(child.age.range)), which(covartemplate$age<min(child.age.range))), which(covartemplate$age<=max(adult.age.range) & covartemplate$age>=min(adult.age.range))),] 
# convert age back to character value so you can use it in a character matrix
covartemplate$age <- as.character(covartemplate$age) 
covartemplate <- as.matrix(covartemplate)
covars <- aaply(covartemplate, 1, function(n){
  # pull out the info you need piece by piece
  sampleID <- n[1] 
  age <- n[3]
  gender <- n[4]
  ethnicity=n[5]
  CELregex <- n[7]
  CELFTP <- n[6]
  study=n[9]
  if(n[10]=="GEO"){
    # for GEO studies, get the CEL file name from the FTP locaiton in the metadata file
    CELfileloc <- file.path(dataDir,study, sub(".gz", "", basename(as.character(CELFTP)))) 
  } else {
    # for ArrayExpress studies, get the CEL file name directly from the CEL file regex in the metadata file
    CELfileloc <- file.path(dataDir, study, CELregex) 
    }
  if(!file.exists(CELfileloc)){
    # check if the CEL file is actually where its supposed to be, if not notify and move on
    print(paste(CELfileloc, "does not exist", sep=" ")) 
    break
    } else {
      # if CEL file is where its supposed to be, add location to covariate table
      return(list(CELfileloc=CELfileloc, ID=sampleID, age=age, gender=gender, ethnicity=ethnicity, study=study)) 
      }
  })
covars <- as.data.frame(covars)
covars$age <- as.numeric(covars$age)
# label samples with their stage as determined by age
# don't need to specify adult ranges, as we limited the dataset above to neonates, children and adults
# was mislabelling 4 year old neonates as adults!
covars$stage <- ifelse(covars$age<=max(neonate.age.range), "neonate", ifelse(covars$age>=min(child.age.range) & covars$age<=max(child.age.range), "child", ifelse(covars$age>=min(adult.age.range) & covars$age<=max(adult.age.range), "adult", NA)))
covars <- covars[order(covars$stage),]
#remove unclassified samples
covars <- covars[!is.na(covars$stage),]


## ----load_data-----------------------------------------------------------
mic.raw <- ReadAffy(filenames=as.character(covars$CELfileloc), phenoData=covars) 
q()

## ----rawQC, eval=FALSE---------------------------------------------------
arrayQualityMetrics(expressionset=mic.raw, outdir=file.path(resultsDir, "QCreport_raw"), force=TRUE, do.logtransform=TRUE, intgroup=c("stage", "study"))


## ----normalize_RMA-------------------------------------------------------
mic.norm.eset <- rma(mic.raw, normalize=TRUE, background=TRUE)


## ----normQC, eval=FALSE--------------------------------------------------
arrayQualityMetrics(expressionset=mic.norm.eset, outdir=file.path(resultsDir, "QCreport_norm"), force=TRUE, do.logtransform=FALSE, intgroup=c("stage", "study"))


## ----drop_outliers, results='asis'---------------------------------------
mic.raw <- mic.raw[,which(!(unlist(pData(mic.raw)$ID) %in% c("GSM844204")))]
mic.norm.eset <- rma(mic.raw, normalize=TRUE, background=TRUE)
save.image(file.path(resultsDir, "RDATA.raw_and_normalized_microarray.data.PBMC.U133Plus2.0" ))


## ----wo.outliers.normQC, eval=FALSE--------------------------------------
arrayQualityMetrics(expressionset=mic.norm.eset, outdir=file.path(resultsDir, "QCreport_norm.wo.outliers"), force=TRUE, do.logtransform=FALSE, intgroup=c("stage", "study"))


## ----print_metadata, results='asis'--------------------------------------
pd <- pData(mic.norm.eset)
pd.gvis <- gvisTable(as.data.frame(apply(pd, 2, as.character)), options=list(width=640))  
print(pd.gvis, "chart")


## ----allexprs------------------------------------------------------------
exprs.all <-exprs(mic.norm.eset)
colnames(exprs.all) <- pData(mic.norm.eset)$ID

probeIDs <- row.names(exprs.all)

symbols <- unlist(mget(rownames(exprs.all), hgu133plus2SYMBOL, ifnotfound=NA))
exprs.all <- cbind(cbind(probeIDs, symbols), exprs.all)
write.table(exprs.all, file.path(resultsDir, "U133_2.0Plus.all.exprs.xls"), sep="\t", row.names=F, col.names=T)


## ----design, results="asis"----------------------------------------------
design <- model.matrix(~ -1+factor(pData(mic.norm.eset)$stage))
# make sure the headings match
colnames(design) <- sub("factor\\(pData\\(mic.norm.eset\\)\\$stage\\)", "", colnames(design))
design.gvis <- gvisTable(as.data.frame(apply(rownames2col(design, "ID"), 2, as.character)), options=list(width=640))  
print(design.gvis, "chart")


## ----contrastmatrix, results='asis'--------------------------------------
contrast.matrix <- makeContrasts(adult-child,adult-neonate,neonate-child, levels=dimnames(design)[[2]])
contrast.gvis <- gvisTable(as.data.frame(apply(rownames2col(contrast.matrix,"contrast"), 2, as.character)), options=list(width=240, height=240))  
print(contrast.gvis, "chart")


## ----linearmodel---------------------------------------------------------
exprs.norm <- exprs(mic.norm.eset)
dimnames(exprs.norm)[[2]] <- as.character(pData(mic.norm.eset)$ID)
fit.exprs <- lmFit(exprs.norm, design) 


## ----contrastfit---------------------------------------------------------
fit2.exprs <- contrasts.fit(fit.exprs, contrast.matrix) 


## ----bayes---------------------------------------------------------------
fit2.exprs <- eBayes(fit2.exprs) 


## ----calcstats, results='hide'-------------------------------------------
stats.exprs <- lapply(seq(1,3,1), function(n) {
  contrast <- dimnames(fit2.exprs$contrasts)$Contrasts[n]
  pd.contrast <- pData(mic.norm.eset)[apply(design[,colnames(design) %in% names(which(abs(contrast.matrix[,contrast])==1))]==1, 1,any),]
  # calculate statistics for this contrast for ALL genes
  stats.contrast <- topTable(fit2.exprs, coef=n, adjust="fdr", sort.by="p", number=nrow(exprs.norm), genelist=row.names(exprs.norm)) 
  # get expression levels for these genes in the samples involved in this contrast
  exprs.contrast <- exprs.norm[,apply(design[,colnames(design) %in% names(which(abs(contrast.matrix[,contrast])==1))]==1, 1,any)]
  # sanity check
  identical(colnames(exprs.contrast), as.vector(unlist(pd.contrast$ID)))
  # rearrange expression values to match order of stats output above
  exprs.contrast <- exprs.contrast[stats.contrast$ID, ]
  # bind columns of stats to columns of sample expressioin values
  stats.exprs.contrast <- cbind(stats.contrast, exprs.contrast)
  # output to file
  write.table(stats.exprs.contrast, file.path(resultsDir,  paste("U133_2.0Plus.stats.exprs", gsub(" ", "", contrast), "xls", sep=".")), sep="\t", row.names=F, col.names=T)
  # sanity check
  identical(colnames(exprs.contrast), as.vector(unlist(pd.contrast$ID)))
  # save for later use in heatmaps
  return(list(contrast=contrast, pd.contrast=pd.contrast, stats.exprs.contrast=stats.exprs.contrast))
  })


## ----heatmaps, fig.width=18, fig.height=18, out.width='100%'-------------
for(n in 1:3){
  # get the names of the stages that are being compared in this comparison
  contrast <- stats.exprs[[n]]$contrast
  pd.contrast <- stats.exprs[[n]]$pd.contrast
  stats.exprs.contrast <- stats.exprs[[n]]$stats.exprs.contrast
  top.stats.exprs.contrast <- subset(stats.exprs.contrast, abs(logFC)>=1 & adj.P.Val<=0.1)
  # remove AFFY control probes
  if (any(grepl("AFFX", top.stats.exprs.contrast$ID))){
    top.stats.exprs.contrast <- top.stats.exprs.contrast[-(grep("AFFX", top.stats.exprs.contrast$ID)),]
    }
  # subset to top 100 (if that many present)
  if (nrow(top.stats.exprs.contrast)>=100){
    top.stats.exprs.contrast <- top.stats.exprs.contrast[1:100,]
    }  
  # extract expression values for samples
  top.exprs.contrast <- top.stats.exprs.contrast[,colnames(top.stats.exprs.contrast) %in% pd.contrast$ID]
  # setup row names for the heatmap, paste probeset ID to gene symbol
  row.names(top.exprs.contrast) <- paste(as.vector(unlist(mget(top.stats.exprs.contrast$ID, hgu133plus2SYMBOL, ifnotfound=NA))), " (" ,top.stats.exprs.contrast$ID, ")", sep="")

  # heatmap annotations
  heatmap.annots <- pd.contrast[,c("ID", "study", "stage", "gender")]
  heatmap.annots <- as.data.frame(apply(heatmap.annots, 2, unlist))
  row.names(heatmap.annots) <- heatmap.annots$ID
  heatmap.annots$ID <- NULL
  # heatmap annotation colors
  study_colors <- c("#FF0000","#00FF00", "#0000FF", cbPalette )
  names(study_colors) <- unique(unlist(pd.contrast$study))
  stage_colors <- c("white", "darkgrey")
  names(stage_colors) <- unique(unlist(pd.contrast$stage))
  gender_colors <- c("cyan", "pink")
  names(gender_colors) <- unique(unlist(pd.contrast$gender))
  ann_colors = list(study = study_colors, stage = stage_colors, gender=gender_colors)
  ## Heatmaps
  # ALL genders 
  pheatmap(as.matrix(top.exprs.contrast), annotation=heatmap.annots, color=rev(brewer.pal(11,"RdBu")), cluster_cols = FALSE,main=paste(contrast,  "- Unclustered", sep=""), show_colnames=F, fontsize=24, fontsize_row=10,annotation_colors=ann_colors)  
  # FEMALE gender  
  top.exprs.contrast.female <- top.exprs.contrast[,which(pd.contrast$gender=="FEMALE")]
  pheatmap(as.matrix(top.exprs.contrast.female), annotation=subset(heatmap.annots,gender=="FEMALE"), cluster_cols = FALSE, color=rev(brewer.pal(11,"RdBu")), main=paste(contrast,"(FEMALE) - Unclustered", sep=" "), show_colnames=F,fontsize=24, fontsize_row=10,annotation_colors=ann_colors)  
  # MALE gender  
  top.exprs.contrast.male <- top.exprs.contrast[,which(pd.contrast$gender=="MALE")]
  pheatmap(as.matrix(top.exprs.contrast.male), annotation=subset(heatmap.annots,gender="MALE"), cluster_cols = FALSE, color=rev(brewer.pal(11,"RdBu")), main=paste(contrast,"(MALE) - Unclustered", sep=" "), show_colnames=F, fontsize=24,fontsize_row=10,annotation_colors=ann_colors)  
  }


## ----heatmaps_all_samples, fig.width=18, fig.height=36-------------------
top.IDs <- unique(unlist(lapply(stats.exprs, function(n) {
  stats.exprs.contrast <- n$stats.exprs.contrast
  top.stats.exprs.contrast <- subset(stats.exprs.contrast, abs(logFC)>=1 & adj.P.Val<=0.1)
  # remove AFFY control probes
  if (any(grepl("AFFX", top.stats.exprs.contrast$ID))){
    top.stats.exprs.contrast <- top.stats.exprs.contrast[-(grep("AFFX", top.stats.exprs.contrast$ID)),]
    }
  # subset to top 100 (if that many present)
  if (nrow(top.stats.exprs.contrast)>=100){
    top.stats.exprs.contrast <- top.stats.exprs.contrast[1:100,]
    } 
  top.stats.exprs.contrast$ID
  })))

# get expression values for these probes
top.exprs.union <- exprs.norm[top.IDs,]

# row labels - add gene symbol to probeset id
row.names(top.exprs.union) <- paste(as.vector(unlist(mget(row.names(top.exprs.union), hgu133plus2SYMBOL, ifnotfound=NA))), " (" ,row.names(top.exprs.union), ")", sep="")
# annotations
pd <- pData(mic.norm.eset)
heatmap.annots <- pd[,c("ID", "study", "stage", "gender")]
heatmap.annots <- as.data.frame(apply(heatmap.annots, 2, unlist))
row.names(heatmap.annots) <- heatmap.annots$ID
heatmap.annots$ID <- NULL
# annotation colors
study_colors <- c("#FF0000","#00FF00", "#0000FF", cbPalette )
names(study_colors) <- unique(unlist(pd$study))
stage_colors <- c("white", "darkgrey", "black")
names(stage_colors) <- unique(unlist(pd$stage))
gender_colors <- c("cyan", "pink")
names(gender_colors) <- unique(unlist(pd$gender))
ann_colors = list(study = study_colors, stage = stage_colors, gender=gender_colors)
## Heatmaps
# Both genders
pheatmap(as.matrix(top.exprs.union), annotation=heatmap.annots, color=rev(brewer.pal(11,"RdBu")), cluster_cols = FALSE, main="All Comparisons, All Samples - Unclustered", show_colnames=F, fontsize=24,fontsize_row=8,annotation_colors=ann_colors)  
# Female gender
top.exprs.union.female <- top.exprs.union[,which(pd$gender=="FEMALE")]
pheatmap(as.matrix(top.exprs.union.female), annotation=subset(heatmap.annots,gender=="FEMALE"), cluster_cols = FALSE, color=rev(brewer.pal(11,"RdBu")), main="All Comparisons, All Female Samples - Unclustered", show_colnames=F,fontsize=24, fontsize_row=8,annotation_colors=ann_colors)  
# Male gender
top.exprs.union.male <- top.exprs.union[,which(pd$gender=="MALE")]
pheatmap(as.matrix(top.exprs.union.male), annotation=subset(heatmap.annots,gender=="MALE"), cluster_cols = FALSE, color=rev(brewer.pal(11,"RdBu")), main="All Comparisons, All male Samples - Unclustered", show_colnames=F, fontsize=24,fontsize_row=8,annotation_colors=ann_colors)   


