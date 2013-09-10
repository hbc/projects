
## @knitr general_libraries
source("http://bioconductor.org/biocLite.R") # BioConductor script necessary for installing new BioC libraries with biocLite()
library(xtable) # table generation for reports
library(plyr) # library for iteratively working with data structures
library(ggplot2) # library for plotting 
library(RColorBrewer) # library for generating color palettes
library(googleVis) # library for presenting tables
source("http://dl.dropboxusercontent.com/u/4253254/Resources/functions.r")


## @knitr general_directories
if (file.exists("/n/hsphS10/hsphfs1/chb/projects/lk_FOY")) {
  baseDir <- "/n/hsphS10/hsphfs1/chb/projects/lk_FOY"
  }  else if (file.exists("/Volumes/home08/jhutchin/consults/lk_FOY/")) {
    baseDir <- "/Volumes/home08/jhutchin/consults/lk_FOY"
    } else {
      baseDir <- "/Volumes/ody/consults/lk_FOY"
      }
dataDir <- file.path(baseDir, "data")
resultsDir <- file.path(baseDir, "results", "WB", "U133Plus2")
metaDir <- file.path(baseDir, "meta", "WB")


## @knitr subset_data
refined.metadata <- read.delim(file.path(metaDir,"unified-metadata-refined.04-17-13.tab"))
U133.Plus2.data <- refined.metadata[which(!is.na(refined.metadata$age) & !is.na(refined.metadata$gender) & !is.na(refined.metadata$CEL_regex) & grepl("GPL570|A-AFFY-44", refined.metadata$platform)),]
write.table(U133.Plus2.data, file.path(metaDir, "unified-metadata-refined_U133Plus2.0.tab"), quote=F, sep="\t", row.names=F, col.names=T)
write.table(U133.Plus2.data, file.path(metaDir, "unified-metadata-refined_U133Plus2.0.xls"), sep="\t", row.names=F, col.names=T)


## @knitr microarray_analysis_libraries
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


## @knitr microarray_analysis_variables
# colorblind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
# age ranges to compare
child.age.range <- c(5,10)
adult.age.range <- c(18,40)
neonate.age.range <- c(0,4)
# THESE CANNOT OVERLAP!


## @knitr covariatedataframe
# U1332.0Plus only for now
covartemplate.file <- "unified-metadata-refined_U133Plus2.0.tab" 
covartemplate <- read.table(file.path(metaDir,covartemplate.file ), header=T, colClasses="character", sep="\t")
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
covars$stage <- ifelse(covars$age<=max(neonate.age.range), "neonate", ifelse(covars$age>=min(child.age.range) & covars$age<=max(child.age.range), "child", ifelse(covars$age>=min(adult.age.range) & covars$age<=max(adult.age.range), "adult", NA)))
covars <- covars[order(covars$stage),]
#remove unclassified samples
covars <- covars[!is.na(covars$stage),]


## @knitr load_data, eval=TRUE
mic.raw <- ReadAffy(filenames=as.character(covars$CELfileloc), phenoData=covars)


## @knitr rawQC, eval=TRUE
#arrayQualityMetrics(expressionset=mic.raw, outdir=file.path(resultsDir, "QCreport_raw"), force=TRUE, do.logtransform=TRUE, intgroup=c("stage", "study"))


## @knitr normalize_RMA, eval=TRUE
mic.norm.eset <- rma(mic.raw, normalize=TRUE, background=TRUE)


## @knitr normQC, eval=TRUE
arrayQualityMetrics(expressionset=mic.norm.eset, outdir=file.path(resultsDir, "QCreport_norm"), force=TRUE, do.logtransform=FALSE, intgroup=c("stage", "study"))


## @knitr drop_outliers, eval=TRUE, results='asis'
mic.raw <- mic.raw[,which(!(unlist(pData(mic.raw)$study) %in% c("GSE18897", "MEXP-2917", "TABM-940", "TABM-666")))]
mic.norm.eset <- rma(mic.raw, normalize=TRUE, background=TRUE)
save.image(file.path(resultsDir, "RDATA.raw_and_normalized_microarray.data_U133Plus2.0" ))


## @knitr normQCwooutliers, eval=TRUE
arrayQualityMetrics(expressionset=mic.norm.eset, outdir=file.path(resultsDir, "QCreport_norm.wo.outliers"), force=TRUE, do.logtransform=FALSE, intgroup=c("stage", "study"))

