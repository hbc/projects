## ----general_libraries---------------------------------------------------
source("http://bioconductor.org/biocLite.R") # BioConductor script necessary for installing new BioC libraries with biocLite()
library(xtable) # table generation for reports
library(plyr) # library for iteratively working with data structures
library(ggplot2) # library for plotting 
library(RColorBrewer) # library for generating color palettes
library(googleVis) # library for presenting tables
library(CHBUtils)

## ----general_directories-------------------------------------------------
if (file.exists("/n/hsphS10/hsphfs1/chb/projects/lk_FOY/continuous_age/normal")) {
  baseDir <- "/n/hsphS10/hsphfs1/chb/projects/lk_FOY/continuous_age/normal"
  }  else if (file.exists("/Volumes/home08/jhutchin/consults/lk_FOY/continuous_age/normal")) {
    baseDir <- "/Volumes/home08/jhutchin/consults/lk_FOY/continuous_age/normal"
    } else {
      baseDir <- "/Volumes/ody/consults/lk_FOY/continuous_age/normal"
      }
dataDir <- file.path(baseDir, "data")
resultsDir <- file.path(baseDir, "results", "PBMC")
metaDir <- file.path(baseDir, "meta", "PBMC")

# colorblind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

## ----subset_data---------------------------------------------------------
refined.metadata <- read.delim(file.path(metaDir,"unified.metadata.refined.tab"))
U133.Plus2.data <- refined.metadata[which(!is.na(refined.metadata$age) & !is.na(refined.metadata$gender) & !is.na(refined.metadata$CEL_regex) & grepl("GPL570|A-AFFY-44", refined.metadata$platform)),]

write.table(U133.Plus2.data, file.path(metaDir, "unified-metadata-refined_U133Plus2.0.tab"), quote=F, sep="\t", row.names=F, col.names=T)
write.table(U133.Plus2.data, file.path(metaDir, "unified-metadata-refined_U133Plus2.0.xls"), sep="\t", row.names=F, col.names=T)

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

## ----covariatedataframe, results='hide'----------------------------------
# U1332.0Plus only for now
covartemplate.file <- "unified-metadata-refined_U133Plus2.0.tab" 
covartemplate <- read.table(file.path(metaDir,covartemplate.file ), header=T, colClasses="character", sep="\t")
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
    print(sampleID)
  if (is.na(n[10])) {
    CELfileloc <- file.path(dataDir, study, CELregex) 
    } else if (n[10]=="GEO"){
      # for GEO studies, get the CEL file name from the FTP locaiton in the metadata file
      CELfileloc <- file.path(dataDir,study, sub(".gz", "", basename(as.character(CELFTP)))) 
      } else {
        # for ArrayExpress or nonGEO studies, get the CEL file name directly from the CEL file regex in the metadata file
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
row.names(covars) <- basename(as.character(covars$CELfileloc))

## ----load_data, eval=FALSE-----------------------------------------------
## mic.raw <- ReadAffy(filenames=as.character(covars$CELfileloc), phenoData=covars)
## save(list="mic.raw", file=file.path(resultsDir, "RDATA.mic.raw"))

## ----rawQC, eval=FALSE---------------------------------------------------
## arrayQualityMetrics(expressionset=mic.raw, outdir=file.path(resultsDir, "QCreport_raw"), force=TRUE, do.logtransform=TRUE, intgroup=c("gender", "study"))

## ----loadesetraw, eval=TRUE, echo=FALSE----------------------------------
# hack to get around reloading raw data into eset
load(file.path(resultsDir, "RDATA.mic.raw"))

## ----exclude_from_raw, eval=TRUE-----------------------------------------
# subset ExpressionSet to studies that are NOT the following
exclude.arrays <- c("GSM340214.CEL","GSM733834.CEL")

mic.raw <- mic.raw[,which(!(row.names(pData(mic.raw)) %in% exclude.arrays))]

## ----normalize_RMA, eval=TRUE--------------------------------------------
mic.norm.eset <- rma(mic.raw,
                     normalize=TRUE,
                     background=TRUE)
save(list="mic.norm.eset", file=file.path(resultsDir, "RDATA.mic.norm.eset"))

## ----normQC, eval=TRUE---------------------------------------------------
arrayQualityMetrics(expressionset=mic.norm.eset, 
                    outdir=file.path(resultsDir, "QCreport_norm"), 
                    force=TRUE,
                    do.logtransform=FALSE, 
                    intgroup=c("gender", "study"))
